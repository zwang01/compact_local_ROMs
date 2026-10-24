function [J, sens, out] = brinkman_stokes_solve_J(par, S, phi)
% BRINKMAN_STOKES_SOLVE_J  Full-order Stokes-Brinkman solve at the design phi
% by one sparse LU factorization, the objective of brinkman_J_eval.m and, when
% requested, the discrete adjoint K' W = dJ/dU (solved with the same factors)
% and the sensitivity dJ/dphi per cell. Used by the tests and the figure
% scripts; the optimization drivers inline the same steps so that they can
% retain the factorization.
[dvec, au, av, dalpha] = brinkman_coeffs(par, S, phi);
K = S.K0 + spdiags(dvec, 0, S.Ntot, S.Ntot);
brk_count_fact;
dK = decomposition(K);
U = dK \ S.F;
[J, Jvisc, dJdU] = brinkman_J_eval(par, S, au, av, U);
Nx = par.Nx; Ny = par.Ny; Nu = S.Nu; Nv = S.Nv;
out = struct('U', U, 'u', reshape(U(1:Nu), Nx+1, Ny), 'v', reshape(U(Nu+1:Nu+Nv), Nx, Ny+1), ...
             'p', reshape(U(Nu+Nv+1:end), Nx, Ny), 'Jvisc', Jvisc);
if nargout < 2, sens = []; return, end
W = dK' \ dJdU;                                  % exact discrete adjoint
sens = brinkman_sens_gather(par, S, dalpha, U, W);
out.W = W;
end
