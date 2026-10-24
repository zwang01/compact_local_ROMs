function F = brk_lu(K, explicit)
% BRK_LU  Sparse LU factorization of the flow saddle-point matrix, returned as
% a struct of solve closures
%     F.solve(x) = K \ x,   F.solvet(x) = K.' \ x.
%
% explicit = false (default): a MATLAB decomposition object, the fastest
%   factorization; used when every system is solved directly (Algorithm 1).
% explicit = true: [L,U,P,Q,R] = lu(K) with P*(R\K)*Q = L*U, whose triangular
%   solves are about three times cheaper than those of decomposition. Used when
%   the factorization is RETAINED as the GMRES preconditioner K_ref of
%   Algorithm 2, where the cost per triangular-solve pair matters.
if nargin < 2, explicit = false; end
if explicit
    [L, U, P, Q, R] = lu(K);
    F.solve  = @(x) Q*(U\(L\(P*(R\x))));
    F.solvet = @(x) R'\(P'*(L'\(U'\(Q'*x))));
else
    dK = decomposition(K);
    F.solve  = @(x) dK \ x;
    F.solvet = @(x) dK' \ x;
end
F.explicit = explicit;
F.n = size(K, 1);
end
