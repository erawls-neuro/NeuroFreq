"""Index source operations across NeuroFreq M-files, App Designer, R and native kernels.

Usage: python tools/audit_arithmetic.py [ROOT] [OUT]. See docs/arithmetic/README.md
for pinned development-only parsers, interpretation and coverage.
"""

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import subprocess
import xml.etree.ElementTree as ET
import zipfile
from bisect import bisect_right
from collections import Counter, defaultdict
from pathlib import Path

import tree_sitter_c
import tree_sitter_cpp
import tree_sitter_r
from miss_hit_core import m_ast as ma
from miss_hit_core.config import Config
from miss_hit_core.errors import Message_Handler
from miss_hit_core.m_lexer import MATLAB_Lexer
from miss_hit_core.m_parser import MATLAB_Parser
from tree_sitter import Language, Parser

C_FAMILIES = {".c", ".h", ".cc", ".cpp", ".cxx", ".hpp", ".hh"}
EXTENSIONS = C_FAMILIES | {".m", ".mlapp", ".r"}
CATEGORIES = {
    "factorization/solve": """chol cholcov cholupdate pinv inv linsolve mldivide mrdivide svd svds
        eig eigs eigval eigvals schur qr lu cond condest rcond rank det null orth lscov lsqminnorm
        mvregress regress robustfit fitlm fitlme fitlmematrix lmer vcov pca cov corr corrcoef
        partialcorr canoncorr crosscov sqrtm expm logm cholesky solve qr.solve""",
    "reduction": """sum cumsum prod cumprod diff mean median geomean harmmean trimmean mad std var norm
        vecnorm max min maxk mink range iqr prctile quantile rms dot cross hypot movsum movmean movmax
        movmin movstd movmedian nanmean nanmedian nanstd nansum skewness kurtosis accumarray trace
        groupmeans varfun splitapply maxabs fmin fmax fabs hypotf hypotl""",
    "elementwise numeric": """plus minus times mtimes rdivide ldivide power mpower uplus uminus abs
        sign sqrt nthroot log log2 log10 log1p exp expm1 pow powf powl logf expf sqrtf floor ceil fix round
        trunc mod rem fmod fmodf real imag angle unwrap conj cos sin tan acos asin atan atan2 cosd
        sind tand acosd asind atan2d sinh cosh tanh asinh acosh atanh erf erfc erfcx erfinv erfcinv
        gamma gammaln gammainc beta betaln betainc factorial nchoosek zscore tiedrank rescale normalize
        deg2rad rad2deg db mag2db db2mag nextpow2 logspace linspace eps double single int8 int16
        int32 int64 uint8 uint16 uint32 uint64 realmin realmax complex cast copysign frexp ldexp
        pmax pmin as.double as.numeric""",
    "numerical algorithm": """fft fft2 fftn ifft ifft2 ifftn fftshift ifftshift conv conv2 convn
        filter filtfilt sosfiltfilt sosfilt fftfilt sinc hilbert cwt wvd spectrogram pwelch tfestimate
        mscohere xcorr interp1 interp2 interpft interpn griddata griddatan interp3 resample detrend
        polyval polyfit gradient trapz integral integral2 quad quadv fminsearch fmincon
        fminbnd fzero optimset lsqcurvefit fit glmfit tcdf tinv normcdf norminv prctile
        binocdf gamcdf betacdf kstest ttest ttest2 butter fir1 firls firpm freqz zp2sos fir2
        bsxfun kron pagemtimes tensorprod mma dotprod interpolate opt shrink rankSVD""",
    "size/axis": """size ndims numel nnz length end height width nrow ncol rows columns ind2sub sub2ind
        reshape permute ipermute squeeze repmat repelem fliplr flipud flip transpose ctranspose
        cat horzcat vertcat diag triu tril blkdiag full sparse eye zeros ones nan inf sparse symamd
        find unique intersect union setdiff ismember sort sortrows seq_len matrix""",
    "numeric predicate": """isfinite isinf isnan isnumeric isfloat isreal issparse isinteger
        any all isequal isequaln isequalwithequalnans isnormal signbit isfinite mxIsNaN
        mxIsInf mxIsFinite ferror isSingular""",
}
CATEGORIES[
    "numeric conversion/axis"
] = """pow2 realpow logspace meshgrid ndgrid pol2cart cart2pol cart2sph sph2cart
    circshift ifftshift shiftdim speye nonzeros typecast bin2dec dec2bin hex2dec dec2hex gather
    gpuArray classUnderlying intmax intmin bitget bitset bitand bitor bitxor bitshift bitcmp logical
    hypotf hypotl sqrtl floorl sinf cosf fabsf fabsl cabs cabsf cpp_int frexp std::hypot"""
CATEGORIES["random sampling"] = (
    "rand randn randi randperm rng randstream datasample randsample normrnd"
)
CATEGORIES[
    "nonnumeric API"
] = """error isempty isfield strcmp struct assert strcmpi isscalar addParameter
    char fullfile exist isstruct ischar fprintf lower sprintf isstring fieldnames isa false cell set strtrim
    fileparts disp get class num2str iscell onCleanup warning cellstr islogical isvector clear validateattributes
    which regexp strjoin isfile javaObject str2double regexprep rethrow rmfield axis title mkdir dir
    strncmpi true parse isfolder addRequired strfind ismatrix path plot3 mfilename addpath datestr cd
    strrep fileread diary orderfields plot xlabel ylabel fclose regexpi ismissing load fopen isrow
    validatestring narginchk startsWith subplot surface strsplit int2str rmdir cla uilabel uibutton
    delete scatter3 close tempname strncmp colorbar patch uieditfield usejava scatter deblank copyfile
    contour table display fwrite getappdata cellfun accumarray arrayfun feval eval evalc str2func"""
# Higher-order calls can invoke arbitrary bodies. Record them explicitly.
CATEGORIES["higher-order API"] = (
    "cellfun arrayfun feval eval evalc str2func bsxfun splitapply"
)
API = {
    name: category for category, names in CATEGORIES.items() for name in names.split()
}
PREDICATES = {"==", "~=", "!=", ">", "<", "<=", ">=", "&&", "||", "&", "|", "~", "!"}
SCALES = {"/", "./", "\\", ".\\", "%", "/=", "%="}
SAFE_CONSTANTS = {
    "eps",
    "pi",
    "realmin",
    "realmax",
    "Inf",
    "NaN",
    "inf",
    "nan",
    "i",
    "j",
}
FIELDS = [
    "source",
    "function",
    "line",
    "column",
    "construct",
    "operation",
    "family",
    "disposition",
    "reason",
    "candidate_policy",
    "source_line",
]


def robust_name(name):
    """Distinguish the project arithmetic module from arbitrary robust estimators."""
    return (
        name.startswith("hnf_rob")
        and not name.startswith("hnf_robust")
        and not name.endswith("_mex")
    )


def classify(file, scope, kind, op, context, literal_divisor=False):
    if kind in {
        "index",
        "range",
        "reshape",
        "constructor",
        "assignment",
        "type/pointer",
        "member",
        "binding",
    }:
        return (
            kind,
            "retain structure",
            "Preserve storage, indexing, sequence or assignment semantics.",
        )
    if kind == "numeric constant":
        return (
            kind,
            "retain primitive",
            "Use the runtime precision/constant as written.",
        )
    if kind in {"call", "handle", "call/index", "macro invocation", "local call"}:
        if kind == "local call":
            return (
                "local dispatch",
                "inspect callee",
                "Local function is present in this source and indexed separately.",
            )
        leaf = re.split(r"::|\.", op)[-1]
        if robust_name(op):
            return (
                "robust helper",
                "delegate",
                "Calls or references the named NeuroFreq implementation; inspect its policy and arguments.",
            )
        if op.endswith("_mex") or kind == "macro invocation":
            return (
                "native kernel",
                "native delegation",
                "Dispatch to a compiled kernel or a recorded local macro body.",
            )
        family = API.get(op, API.get(leaf))
        if op.startswith(("nf_", "hnf_")) and not family:
            return (
                "NeuroFreq dispatch",
                "inspect callee",
                "Arithmetic in the target is recorded at its definition.",
            )
        if family in {"nonnumeric API", "higher-order API", "random sampling"}:
            return (
                family,
                "retain interface",
                "Preserve runtime dispatch, parameter checks, random stream or display/I/O contract.",
            )
        if family == "size/axis" or family == "numeric predicate":
            return (
                family,
                "retain structure",
                "Preserve the requested axes, ordering or exact predicate.",
            )
        if family in {"factorization/solve", "numerical algorithm"}:
            return (
                family,
                "retain algorithm",
                "Library or method owns an estimator, decomposition, transform or optimization; changing tolerances/regularization changes its contract.",
            )
        if family in {"reduction", "elementwise numeric", "numeric conversion/axis"}:
            if "/robustarithmetic/" in file or robust_name(scope):
                return (
                    family,
                    "helper primitive",
                    "Part of the named implementation; re-entering the wrapper would change semantics or recurse.",
                )
            return (
                family,
                "retain numerical primitive",
                "Preserve the native domain, NaN policy, reduction dimension and type; inspect the source expression for scale or singularity assumptions.",
            )
        if op.startswith(("mx", "mex")) or kind == "handle":
            return (
                "API/handle",
                "retain interface",
                "Runtime API or callback; no scalar replacement established.",
            )
        return (
            "other call/index",
            "unresolved at static analysis",
            "Includes plotting, I/O, local/vendor functions, dynamic handles and unresolved indexing. The source line locates the use.",
        )
    if op in PREDICATES or kind == "predicate":
        return (
            "predicate/logical",
            "retain predicate",
            "The comparison or logical result controls inclusion, branching, masks or validation.",
        )
    if kind in {"transpose", "pointer"}:
        return (
            kind,
            "retain primitive",
            "Preserve conjugation, shape or address semantics.",
        )
    if context == "index":
        return (
            "index arithmetic",
            "retain index arithmetic",
            "Operand appears inside an array subscript; rounding or clamping can select a different element.",
        )
    if "/robustarithmetic/" in file or robust_name(scope):
        return (
            "arithmetic",
            "helper primitive",
            "This is the implementation of the operation or its support; review its documented domain and failure policy.",
        )
    if file.startswith("tests/"):
        return (
            "arithmetic",
            "test/reference",
            "Independent test data, expected results or tolerance; wrapping can hide the regression under test.",
        )
    if op in SCALES:
        if literal_divisor:
            return (
                "division/solve",
                "retain scale",
                "Division by a nonzero numeric literal; a denominator floor would add no protection here.",
            )
        return (
            "division/solve",
            "retain; review domain",
            "A fixed denominator floor or pseudoinverse could change units, mask missingness or change a model rank. Inspect domain/conditioning in context.",
        )
    if op in {"^", ".^", "**"}:
        return (
            "power",
            "retain; review domain",
            "Preserve the scalar/elementwise/matrix exponent contract and complex domain. Inspect extreme magnitudes at this site.",
        )
    if op in {"*", ".*", "*="}:
        return (
            "multiplication",
            "retain formula",
            "Operator may be scaling or a product; wrapping inputs or saturating the result can change the estimator.",
        )
    if op in {"+", "-", "+=", "-=", "++", "--"}:
        return (
            "addition/subtraction",
            "retain formula",
            "Preserve units, exact index/count steps or signed residuals. Cancellation depends on operand scale.",
        )
    return (
        "other operator",
        "retain primitive",
        "Bit, logical, assignment variant or language-specific operator; inspect source context.",
    )


# Candidate policies are audit prompts, not a source-to-source rewrite rule.
CANDIDATES = {
    "solve": "hnf_robsolve for matrix systems; verify rank, units, residual and regularization.",
    "division": "hnf_robdiv only after specifying behavior at zero, NaN and an unresolved denominator.",
    "sum": "hnf_robsum / hnf_robmean when cancellation or overflow is plausible; preserve dimension and missingness.",
    "moments": "hnf_robvar, hnf_robcov, hnf_robcorr, or hnf_robzscore with the existing sample weights and degrees of freedom.",
    "power": "hnf_robabs2 / hnf_robpower for magnitude squared; general matrix or signed powers need their original definition.",
    "domain": "Check domain and representable range; preserve complex inputs and deliberate infinities.",
}


def candidate(family, op, disposition):
    if disposition not in {
        "retain; review domain",
        "retain algorithm",
        "retain numerical primitive",
    }:
        return ""
    if family == "division/solve":
        return (
            CANDIDATES["solve"] + " " + CANDIDATES["division"]
            if op in {"/", "\\"}
            else CANDIDATES["division"]
        )
    if family == "power":
        return CANDIDATES["power"]
    if op in {"chol", "cholcov"}:
        return "hnf_robchol only with explicit jitter policy and effective matrix reporting."
    if op in {"inv", "pinv", "mldivide", "mrdivide", "linsolve"}:
        return "hnf_robpinv or " + CANDIDATES["solve"]
    if op in {"sum", "mean", "prod", "cumsum", "nanmean", "nansum"}:
        return (
            CANDIDATES["sum"]
            if op not in {"prod", "cumsum"}
            else "Keep cumulative/product contract; inspect accumulation range."
        )
    if op in {"var", "std", "zscore", "cov", "corr", "corrcoef"}:
        return CANDIDATES["moments"]
    if op in {"exp", "log", "log10", "log2", "sqrt", "atanh", "acos", "asin"}:
        return CANDIDATES["domain"]
    if op in {"norm", "rms", "dot"}:
        return "Use a scaled norm/RMS or compensated dot product when the native calculation overflows or cancels."
    return ""


class Source:
    def __init__(self, file, contents):
        self.file = file
        self.lines = contents.splitlines()
        self.rows = []

    def add(self, scope, line, col, kind, op, context="", literal_divisor=False):
        summary = classify(self.file, scope, kind, op, context, literal_divisor)
        self.rows.append(
            [
                self.file,
                scope,
                line,
                col,
                kind,
                op,
                *summary,
                candidate(summary[0], op, summary[1]),
                self.lines[line - 1].strip(),
            ]
        )


def symbol(node):
    # str(Name) is correct for simple and dotted names. Compound subscripts
    # can contain arbitrary expression nodes whose default repr embeds a pointer.
    return re.sub(r"<[^<>]* object at 0x[\da-fA-F]+>", "<expression>", str(node))


class MatlabVisitor(ma.AST_Visitor):
    def __init__(self, source, collect_bindings):
        self.source = source
        self.collect_bindings = collect_bindings
        self.scopes = ["<script>"]
        self.bound = defaultdict(set)
        self.functions = set()

    def root_name(self, node):
        name = re.match(r"[A-Za-z]\w*", symbol(node))
        return name.group() if name else ""

    def is_bound(self, node):
        return any(self.root_name(node) in self.bound[scope] for scope in self.scopes)

    def index_context(self, node):
        child = node
        parent = node.n_parent
        while parent is not None and isinstance(parent, ma.Expression):
            if (
                isinstance(parent, (ma.Reference, ma.Cell_Reference))
                and child in parent.l_args
                and (
                    isinstance(parent, ma.Cell_Reference)
                    or self.is_bound(parent.n_ident)
                )
            ):
                return "index"
            child, parent = parent, parent.n_parent
        return ""

    def visit(self, node, parent, relation):
        if isinstance(node, ma.Function_Definition):
            self.scopes.append(symbol(node.n_sig.n_name))
        scope = self.scopes[-1]
        if self.collect_bindings:
            if isinstance(node, ma.Function_Definition):
                sig = node.n_sig
                self.functions.add(symbol(sig.n_name))
                targets = (sig.l_inputs or []) + (sig.l_outputs or [])
            elif isinstance(node, ma.Simple_Assignment_Statement):
                targets = [node.n_lhs]
            elif isinstance(node, ma.Compound_Assignment_Statement):
                targets = node.l_lhs or []
            elif isinstance(node, ma.Lambda_Function):
                targets = node.l_parameters
            elif isinstance(node, ma.For_Loop_Statement):
                targets = [node.n_ident]
            else:
                return
            self.bound[scope].update(self.root_name(x) for x in targets)
            return

        literal_divisor = False
        if isinstance(node, (ma.Unary_Operation, ma.Binary_Operation)):
            op = node.t_op.value
            kind = "unary" if isinstance(node, ma.Unary_Operation) else "binary"
            if op in {"'", ".'"}:
                kind = "transpose"
            if isinstance(node, ma.Binary_Operation) and isinstance(
                node.n_rhs, ma.Number_Literal
            ):
                try:
                    literal_divisor = float(str(node.n_rhs)) != 0
                except ValueError:
                    pass
        elif isinstance(node, ma.Identifier) and str(node) == "~":
            kind, op = "binding", "~"
        elif isinstance(node, ma.Range_Expression):
            for token in (node.t_first_colon, node.t_second_colon):
                if token:
                    loc = token.location
                    self.source.add(scope, loc.line, loc.col_start + 1, "range", ":")
            return
        elif isinstance(node, ma.Reshape):
            kind, op = "reshape", ":"
        elif isinstance(
            node, (ma.Simple_Assignment_Statement, ma.Compound_Assignment_Statement)
        ):
            kind, op = "assignment", "="
        elif isinstance(node, (ma.Matrix_Expression, ma.Cell_Expression)):
            kind, op = (
                "constructor",
                "[]" if isinstance(node, ma.Matrix_Expression) else "{}",
            )
        elif isinstance(node, ma.Function_Pointer):
            kind, op = "handle", symbol(node.n_name)
        elif isinstance(node, ma.Function_Call):
            kind, op = "call", symbol(node.n_name)
        elif isinstance(node, (ma.Reference, ma.Cell_Reference)):
            op = symbol(node.n_ident)
            kind = "index" if isinstance(node, ma.Cell_Reference) else "call/index"
            if isinstance(node, ma.Reference):
                if self.is_bound(node.n_ident):
                    kind = "index"
                elif op in self.functions:
                    kind = "local call"
                elif op in API or op.startswith(("nf_", "hnf_")):
                    kind = "call"
        elif (
            isinstance(node, ma.Identifier)
            and str(node) in SAFE_CONSTANTS
            and not self.is_bound(node)
        ):
            if not isinstance(
                parent, (ma.Reference, ma.Selection, ma.Function_Pointer)
            ):
                kind, op = "numeric constant", str(node)
            else:
                return
        else:
            return
        loc = node.loc()
        self.source.add(
            scope,
            loc.line,
            loc.col_start + 1,
            kind,
            op,
            self.index_context(node),
            literal_divisor,
        )

    def visit_end(self, node, parent, relation):
        if isinstance(node, ma.Function_Definition):
            self.scopes.pop()


def matlab(source, text):
    config = Config()
    messages = Message_Handler("lint")
    messages.register_file(source.file)
    tree = MATLAB_Parser(
        messages, MATLAB_Lexer(config.language, messages, text, source.file), config
    ).parse_file()
    visitor = MatlabVisitor(source, True)
    tree.visit(None, visitor, "Root")
    visitor.collect_bindings = False
    tree.visit(None, visitor, "Root")
    # Re-lex independently of AST traversal. Every arithmetic/logical operator
    # and colon in source, including ignored "~" argument placeholders, must
    # appear at its exact token position. Comments and strings are not tokens.
    records = {(row[2], row[3], row[5]) for row in source.rows}
    lexer = MATLAB_Lexer(config.language, messages, text, source.file)
    while (token := lexer.token()) is not None:
        if token.kind in {"OPERATOR", "COLON"}:
            op = ":" if token.kind == "COLON" else token.value
            key = (token.location.line, token.location.col_start + 1, op)
            if key not in records:
                raise RuntimeError(f"{source.file}: unaccounted source token {key}")
    if messages.errors:
        raise RuntimeError(f"{messages.errors} MATLAB parser diagnostics")


def node_text(node):
    return node.text.decode("utf-8")


def function_name(node):
    # Unwrap nested C declarators (pointers, function declarators, attributes).
    while node is not None:
        if node.type in {"identifier", "field_identifier", "qualified_identifier"}:
            return node_text(node)
        node = node.child_by_field_name("declarator")
    return "<anonymous>"


def walk_tree(parser, source, text, language):
    tree = parser.parse(text)
    macros = {
        node_text(n.child_by_field_name("name"))
        for n in tree.root_node.children
        if n.type == "preproc_function_def"
    }
    starts = [0] + [match.end() for match in re.finditer(b"\n", text)]

    def walk(node, scope, subscript=False, offset=0):
        if node.type == "ERROR" or node.is_missing:
            raise RuntimeError(
                f"{source.file}:{node.start_point} {node.type}: {node.text[:80]!r}"
            )
        kind = node.type

        def field(name):
            return node.child_by_field_name(name)

        def record(position, construct, op, constant=False):
            # Tree-sitter columns count bytes. Map the byte offset back to a
            # character column, including macro bodies parsed out of context.
            absolute = position.start_byte + offset
            line = bisect_right(starts, absolute)
            start = starts[line - 1]
            column = len(text[start:absolute].decode("utf-8")) + 1
            source.add(
                scope,
                line,
                column,
                construct,
                op,
                "index" if subscript else "",
                constant,
            )

        if kind == "function_definition":
            if language == "R":
                scope = (
                    node_text(node.parent.child_by_field_name("lhs"))
                    if node.parent.type == "binary_operator"
                    else "<anonymous>"
                )
            else:
                scope = function_name(field("declarator"))
        if kind in {"preproc_function_def", "preproc_def"}:
            body = field("value")
            if body:
                # Parse the literal replacement text once, including line
                # continuations; do not multiply operations by macro call count.
                prefix = b"void audit_macro(void){\n"
                payload = body.text.replace(b"\\\r\n", b"  \r\n").replace(
                    b"\\\n", b" \n"
                )
                unit = parser.parse(prefix + payload + b";\n}")
                wrapped = unit.root_node.named_children[0].child_by_field_name("body")
                for part in wrapped.named_children:
                    walk(
                        part,
                        node_text(field("name")),
                        False,
                        body.start_byte - len(prefix),
                    )
            return
        if kind in {
            "binary_expression",
            "unary_expression",
            "update_expression",
            "assignment_expression",
            "pointer_expression",
            "binary_operator",
            "unary_operator",
            "namespace_operator",
            "extract_operator",
        }:
            operator = field("operator")
            if operator is None:
                operator = next(c for c in node.children if not c.is_named)
            op = node_text(operator)
            construct = (
                "unary"
                if kind in {"unary_expression", "unary_operator", "update_expression"}
                else "binary"
            )
            if kind == "pointer_expression":
                construct = "pointer"
            if kind in {"namespace_operator", "extract_operator"}:
                construct = "member"
            if op in {"=", "<-", "->", "<<-", "->>"}:
                construct = "assignment"
            rhs = field("right") or field("rhs")
            constant = False
            if rhs and rhs.type in {"number_literal", "float", "integer"}:
                try:
                    constant = float(node_text(rhs).rstrip("fFlLuU")) != 0
                except ValueError:
                    pass
            record(operator, construct, op, constant)
        elif kind in {"call_expression", "call"}:
            ident = field("function")
            callee = node_text(ident)
            record(ident, "macro invocation" if callee in macros else "call", callee)
        elif kind in {"subscript_expression", "subset", "subset2"}:
            token = next((c for c in node.children if c.type in {"[", "[["}), node)
            record(token, "index", node_text(token) if token != node else "[]")
        elif kind == "cast_expression":
            record(node, "type/pointer", node_text(field("type")))
        elif kind == "conditional_expression":
            record(next(c for c in node.children if c.type == "?"), "predicate", "?:")
        for child in node.children:
            if child.is_named or child.is_missing:
                inner = (
                    subscript
                    or (kind == "subscript_expression" and child == field("index"))
                    or (kind in {"subset", "subset2"} and child == field("arguments"))
                )
                walk(child, scope, inner, offset)

    walk(tree.root_node, "<file>")


def source_paths(root):
    try:
        files = subprocess.check_output(
            ["git", "-C", str(root), "ls-files", "-z", "-co", "--exclude-standard"],
            stderr=subprocess.DEVNULL,
        )
        paths = (root / item.decode() for item in files.split(b"\0") if item)
    except (FileNotFoundError, subprocess.CalledProcessError):
        paths = (
            p
            for p in root.rglob("*")
            if not {".git", "_nf_local", "__pycache__"}.intersection(p.parts)
        )
    return sorted({p for p in paths if p.suffix.lower() in EXTENSIONS and p.is_file()})


def csv_file(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument(
        "root", nargs="?", type=Path, default=Path(__file__).resolve().parents[1]
    )
    cli.add_argument("out", nargs="?", type=Path)
    args = cli.parse_args()
    root = args.root.resolve()
    out = args.out or root / "docs" / "arithmetic"
    out.mkdir(parents=True, exist_ok=True)
    parsers = {
        "C": Parser(Language(tree_sitter_c.language())),
        "C++": Parser(Language(tree_sitter_cpp.language())),
        "R": Parser(Language(tree_sitter_r.language())),
    }
    coverage = []
    arithmetic = Counter()
    helper_sites = Counter()
    counts = Counter()
    treehash = hashlib.sha256()
    with (
        (out / "operations.csv.gz").open("wb") as raw,
        gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz,
        io.TextIOWrapper(gz, encoding="utf-8", newline="") as stream,
    ):
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(FIELDS)
        for path in source_paths(root):
            name = path.relative_to(root).as_posix()
            original = path.read_bytes()
            digest = hashlib.sha256(original).hexdigest()
            treehash.update((name + "\0" + digest + "\n").encode())
            suffix = path.suffix.lower()
            if suffix == ".mlapp":
                with zipfile.ZipFile(io.BytesIO(original)) as archive:
                    xml = ET.fromstring(archive.read("matlab/document.xml"))
                text = "".join(xml.itertext()).strip()
                name += "!matlab/document.xml"
                language = "MATLAB app"
            else:
                text = original.decode("utf-8-sig")
                language = (
                    "MATLAB"
                    if suffix == ".m"
                    else "R"
                    if suffix == ".r"
                    else "C"
                    if suffix in {".c", ".h"}
                    else "C++"
                )
            source = Source(name, text)
            if language.startswith("MATLAB"):
                matlab(source, text.replace("\r\n", "\n").replace("\r", "\n"))
            else:
                walk_tree(parsers[language], source, text.encode("utf-8"), language)
            source.rows.sort(key=lambda row: (row[2], row[3], row[4], row[5]))
            for row in source.rows:
                ln, col, role, op = row[2:6]
                if role in {
                    "unary",
                    "binary",
                    "range",
                    "reshape",
                    "transpose",
                    "pointer",
                    "predicate",
                    "assignment",
                    "binding",
                    "member",
                }:
                    expected = "?" if op == "?:" else op
                    if not source.lines[ln - 1][col - 1 :].startswith(expected):
                        raise RuntimeError(f"{name}:{ln}:{col} does not contain {op!r}")
                if role in {"unary", "binary"} and op in {
                    "+",
                    "-",
                    "*",
                    ".*",
                    "/",
                    "./",
                    "\\",
                    ".\\",
                    "^",
                    ".^",
                    "%",
                    "++",
                    "--",
                    "+=",
                    "-=",
                    "*=",
                    "/=",
                    "%=",
                    "**",
                }:
                    arithmetic[op] += 1
                if row[7] == "delegate" and not name.startswith(
                    ("tests/", "helpers/robustarithmetic/")
                ):
                    helper_sites[(op, role)] += 1
            writer.writerows(source.rows)
            counts.update(
                (
                    language,
                    r[4],
                    r[5] if r[4] not in {"index", "call/index"} else "(all)",
                    r[6],
                    r[7],
                )
                for r in source.rows
            )
            coverage.append(
                [name, language, digest, len(source.lines), len(source.rows)]
            )
    csv_file(
        out / "coverage.csv",
        ["source", "language", "sha256", "lines", "rows"],
        coverage,
    )
    csv_file(
        out / "operation_counts.csv",
        ["language", "construct", "operation", "family", "disposition", "sites"],
        [
            (*key, n)
            for key, n in sorted(counts.items(), key=lambda pair: (-pair[1], pair[0]))
        ],
    )
    print(
        f"{len(coverage)} files, {sum(counts.values())} records, source SHA-256 {treehash.hexdigest()}"
    )
    snapshot = {
        "source_tree_sha256": treehash.hexdigest(),
        "files": len(coverage),
        "records": sum(counts.values()),
        "by_language": {
            lang: {
                "files": sum(c[1] == lang for c in coverage),
                "records": sum(c[4] for c in coverage if c[1] == lang),
            }
            for lang in sorted({c[1] for c in coverage})
        },
    }
    snapshot["arithmetic_operators"] = dict(arithmetic.most_common())
    snapshot["runtime_arithmetic_helper_sites"] = {
        name + " (" + role + ")": n for (name, role), n in helper_sites.most_common()
    }
    summary = [
        "# Snapshot counts",
        "",
        f"Scanned **{len(coverage)} source files** and recorded **{sum(counts.values()):,} source sites**. "
        f"Among them, **{sum(arithmetic.values()):,}** are explicit arithmetic operators. Other sites include numerical calls, "
        "transposes, masks, indexing, assignments, and dispatch. Counts measure source occurrences; a loop body contributes once.",
        "",
        "| Language | Files | Source sites |",
        "| --- | ---: | ---: |",
    ]
    for language, totals in snapshot["by_language"].items():
        summary.append(f"| {language} | {totals['files']} | {totals['records']:,} |")
    summary.extend(
        [
            "",
            "## Explicit arithmetic operators",
            "",
            "Counts include tests and index calculations. In MATLAB, `*` and `/` may have scalar or matrix inputs. A dot marks elementwise arithmetic.",
            "",
            "| Operator | Sites |",
            "| --- | ---: |",
        ]
    )
    summary.extend(f"| `{op}` | {n:,} |" for op, n in arithmetic.most_common())
    summary.extend(
        [
            "",
            "## Calls to robust arithmetic from outside the arithmetic module and tests",
            "",
            "| Helper or handle | Sites |",
            "| --- | ---: |",
        ]
    )
    summary.extend(
        f"| `{op}` ({role}) | {n:,} |" for (op, role), n in helper_sites.most_common()
    )
    summary.extend(
        [
            "",
            f"Source tree SHA-256 `{treehash.hexdigest()}`. Individual input digests are in `coverage.csv`.",
            "",
        ]
    )
    (out / "snapshot.md").write_text("\n".join(summary), encoding="utf-8")
    (out / "snapshot.json").write_text(
        json.dumps(snapshot, indent=2) + "\n", encoding="utf8"
    )


if __name__ == "__main__":
    main()
