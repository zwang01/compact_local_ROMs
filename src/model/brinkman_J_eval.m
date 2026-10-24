function [J, Jvisc, dJdU, Jdrag] = brinkman_J_eval(par, S, au, av, U)
% BRINKMAN_J_EVAL  Objective and its exact gradient with respect to the
% discrete state U (full-order or reduced). par.objective selects
%   'dissipation' (default): J = int nu|grad u|^2 + alpha|u|^2   [eq:Jflow],
%   'drag'                 : J = int alpha u_x                   [eq:dragdef].
% Jvisc (viscous part of the dissipation) and Jdrag are returned for reporting.
h = par.h; nu = par.nu;
Nu = S.Nu; Nv = S.Nv;
uv = U(1:Nu+Nv);
u  = U(1:Nu); v = U(Nu+1:Nu+Nv);
g  = S.E*uv + S.e0;
Jvisc = h^2 * nu*(g.'*g);
Jdiss = Jvisc + h^2 * ( sum(au(:).*u.^2) + sum(av(:).*v.^2) );
Jdrag = h^2 * sum(au(:).*u);          % x-force on the solid: the Brinkman reaction

if isfield(par,'objective') && strcmp(par.objective,'drag')
    J = Jdrag;
    if nargout > 2
        dJdU = [h^2*au(:); zeros(S.Nv,1); zeros(S.Np,1)];
    end
else
    J = Jdiss;
    if nargout > 2
        dJduv = 2*h^2 * ( nu*(S.E.'*g) + [au(:).*u; av(:).*v] );
        dJdU  = [dJduv; zeros(S.Np,1)];
    end
end
end
