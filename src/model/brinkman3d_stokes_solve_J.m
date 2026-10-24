function [J, sens, out] = brinkman3d_stokes_solve_J(par, S, phi)
% BRINKMAN3D_STOKES_SOLVE_J  3D analogue of brinkman_stokes_solve_J.m: one
% sparse LU solve of the Stokes-Brinkman system at the design phi, the
% objective of brinkman3d_J_eval.m and, when requested, the exact discrete
% adjoint and the sensitivity dJ/dphi. Used by the tests and the figures.
[dvec, au, av, aw, dalpha] = brinkman3d_coeffs(par, S, phi);
K = S.K0 + spdiags(dvec, 0, S.Ntot, S.Ntot);
brk_count_fact;
dK = decomposition(K);
U = dK \ S.F;
[J, Jvisc, dJdU, Jdrag] = brinkman3d_J_eval(par, S, au, av, aw, U);
Nx = par.Nx; Ny = par.Ny; Nz = par.Nz; Un = S.Un; Vn = S.Vn; Wn = S.Wn;
out = struct('U', U, 'u', reshape(U(1:Un),Nx+1,Ny,Nz), 'v', reshape(U(Un+1:Un+Vn),Nx,Ny+1,Nz), ...
             'w', reshape(U(Un+Vn+1:Un+Vn+Wn),Nx,Ny,Nz+1), 'p', reshape(U(Un+Vn+Wn+1:end),Nx,Ny,Nz), ...
             'Jvisc', Jvisc, 'Jdrag', Jdrag);
if nargout < 2, sens = []; return, end
W = dK' \ dJdU;
sens = brinkman3d_sens_gather(par, S, dalpha, U, W);
out.W = W;
end
