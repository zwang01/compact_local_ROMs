function [J, Jvisc, dJdU, Jdrag] = brinkman3d_J_eval(par, S, au, av, aw, U)
% BRINKMAN3D_J_EVAL  3D analogue of brinkman_J_eval.m: objective ('drag',
% int alpha u_x, or 'dissipation') and its exact gradient with respect to U.
h = par.h; nu = par.nu; Un = S.Un; Vn = S.Vn; Wn = S.Wn;
uvw = U(1:Un+Vn+Wn);
u = U(1:Un); v = U(Un+1:Un+Vn); w = U(Un+Vn+1:Un+Vn+Wn);
g = S.E*uvw;
Jvisc = h^3*nu*(g.'*g);
Jdiss = Jvisc + h^3*(sum(au(:).*u.^2) + sum(av(:).*v.^2) + sum(aw(:).*w.^2));
Jdrag = h^3*sum(au(:).*u);
if isfield(par,'objective') && strcmp(par.objective,'drag')
    J = Jdrag;
    if nargout > 2, dJdU = [h^3*au(:); zeros(Vn+Wn+S.Pn,1)]; end
else
    J = Jdiss;
    if nargout > 2
        dJdU = [2*h^3*(nu*(S.E.'*g) + [au(:).*u; av(:).*v; aw(:).*w]); zeros(S.Pn,1)];
    end
end
end
