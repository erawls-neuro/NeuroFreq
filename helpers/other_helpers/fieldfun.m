function S = fieldfun(fun,S1,Sn)
%FIELDFUN Apply a function to the matching fields of structures.
%   S = FIELDFUN(FUN,S1,...,SN) passes the values for each field in
%   structures S1,...,SN to function FUN and returns the result in scalar
%   structure S. Each field F in S, i.e. S.(F), is the result of 
%   FUN(S1.(F),...,SN.(F)).
% 
%   - FUN is either a function handle or a text scalar giving the 
%       function's name or operator. The function must take N input 
%       arguments, i.e. equal to the number of input structures, and 
%       return a value. FUN is called once for each field F in S1,...,SN.
%   - S1 is a scalar structure or array of structures.
%   - S2,...,SN are scalar structures, arrays of structures, or variables 
%       of any other class. If S2,...,SN are structures, their fields and 
%       field order must match S1. If S2,...,SN are not structures, then 
%       they are converted to structures, with the value being assigned to 
%       each field.
%   - S is a scalar structure with the same fields and field order as S1.
%       The value of each field F in S, i.e. S.(F), is the result of 
%       FUN(S1.(F),...,SN.(F)).
% 
%   Examples:
%       S(1) = struct( 'name', "John", 'age', 30, 'vaccinated', true );
%       S(2) = struct( 'name', "Jane", 'age', 80, 'vaccinated', false );
%       FUN = @horzcat;
%       fieldfun( FUN, S )
%       % Returns struct( 'name', ["John" "Jane"], 'age', [30 80], ...
%       %     'vaccinated', [true false] ).
%       FUN = @(varargin) replace( string( varargin ), ...
%           alphanumericsPattern(1), "#" );
%       fieldfun( FUN, S )
%       % Returns struct( 'name', ["####","####"], 'age', ["##","##"], ...
%       %     'vaccinated', ["####","#####"] ).
%         
%       IsWorkDone = struct( 'methods', true, 'results', true );
%       IsReportWritten = struct( 'methods', true, 'results', true );
%       IsChecked = struct( 'methods', true, 'results', false );
%       FUN = @(varargin) all([varargin{:}]);
%       fieldfun( FUN, IsWorkDone, IsReportWritten, IsChecked )
%       % Returns struct( 'methods', true, 'results', false )
%         
%       IsRaining = struct( 'Sun', false, 'Mon', false, 'Tue', true );
%       isHoliday = true;
%       fieldfun( "&", fieldfun( "~", IsRaining ), isHoliday )
%       % Returns struct( 'Sun', true, 'Mon', true, 'Tue', false )
% 
%       For more detailed examples see examples.mlx/examples.pdf.
%    
%   Created in 2022b. Compatible with 2019b and later. Compatible with all 
%   platforms. Please cite George Abrahams 
%   https://github.com/WD40andTape/fieldfun.
% 
%   See also STRUCTFUN, CELLFUN, ARRAYFUN, VARFUN, ROWFUN, FUNCTION_HANDLE

%   Published under MIT License (see LICENSE.txt).
%   Copyright (c) 2022 George Abrahams.
%   - https://github.com/WD40andTape/
%   - https://www.linkedin.com/in/georgeabrahams/

    arguments
        fun {mustBeFunctionHandleOrTextScalar}
        S1 {mustBeA(S1,'struct')}
    end
    arguments (Repeating)
        Sn {mustHaveSameFieldsIfStruct(Sn,S1)}
    end

    if istextscalar(fun)
        fun = str2func(string(fun));
    end
    fields = fieldnames(S1);

    isConvertToStruct = ~cellfun(@isstruct,Sn);
    if any(isConvertToStruct)
        Sn(isConvertToStruct) = num2cell( cell2struct( ...
            repmat(Sn(isConvertToStruct),numel(fields),1), fields, 1 ) );
    end

    S = [S1 Sn{:}];
    S = cell2struct( cellfun( @(x)fun(S.(x)), fields, ...
        'UniformOutput', false), fields, 1 );
end

function mustBeFunctionHandleOrTextScalar(a)
    if ~(isa(a,'function_handle') || istext(a))
        eidType = 'fieldfun:Validators:NotFunctionHandleOrText';
        msgType = 'Value must be a function handle or function name.';
        throwAsCaller(MException(eidType,msgType))
    end
    if ~(isscalar(a) || istextscalar(text))
        eidType = 'fieldfun:Validators:NotScalarOrTextScalar';
        msgType = 'Value must be a scalar or text scalar.';
        throwAsCaller(MException(eidType,msgType))
    end
end
function mustHaveSameFieldsIfStruct(A,B)
    if isstruct(A) && ~isequal(fieldnames(A),fieldnames(B))
        eidType = 'fieldfun:Validators:InconsistentFieldNames';
        msgType = ['Scalar structure inputs must have the same fields,' ...
            ' in the same order.'];
        throwAsCaller(MException(eidType,msgType))
    end
end

function tf = istext(text)
    tf = ischar(text) | iscellstr(text) | isstring(text);
end
function tf = istextscalar(text)
    tf = ( isstring(text) & isscalar(text) ) | ...
        ( ischar(text) & ( isrow(text) | isequal(size(text),[0 0]) ) );
end