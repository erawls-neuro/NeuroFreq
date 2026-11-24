function [H, scale, defect] = hnf_la_hermitian(A, symmetry, symmetryTolerance)
% Internal: scale before Hermitian averaging to prevent A+A' overflow.
hnf_la_check(A);
if size(A,1) ~= size(A,2), error('hnf:Shape','Expected a square matrix.'); end
if ~any(strcmp(symmetry,{'check','symmetrize'}))
    error('hnf:InvalidOption','symmetry must be check or symmetrize.');
end
hnf_la_scalar(symmetryTolerance,'symmetryTolerance',0,1);
scale = hnf_la_scale(A);
H = A./scale;
defect = norm(H-H','fro')/max(norm(H,'fro'),realmin(class(A)));
if strcmp(symmetry,'check') && defect > symmetryTolerance
    error('hnf:NotHermitian','Matrix is not Hermitian within symmetryTolerance.');
end
H = H/2 + H'/2;
end
