function [phi, info] = ac_update_vol(phi, sens, freeze, phiFrozen, Vcells, ac)
% AC_UPDATE_VOL  One design update of the paper (Algorithm 1, lines 7-10):
% ac.Kin inner steps of the gated, volume-constrained Allen-Cahn flow
% [eq:ac_update_vol], each solving the semi-implicit system
%   (I - tau*hat_eps^2*L_h) phi_new = r0 + lambda*r1,
%   r0 = phi + tau*(-F'(phi) - beta*G'(phi).*ghat),   r1 = tau*G'(phi),
% with the double well F(phi) = phi^2(1-phi)^2/4, F'(phi) = phi(1-phi)(1/2-phi),
% the gate G(phi) = phi^3(6phi^2-15phi+10)/30, G'(phi) = phi^2(1-phi)^2, and
% ghat = sens/||sens|| the normalized sensitivity, held fixed over the inner
% steps. L_h is the cell-centred Laplacian and hat_eps = eps/(3 sqrt 2). The update is affine
% in lambda, phi_new = phi_a + lambda*phi_b, with phi_a, phi_b from two solves
% with the prefactorized matrix ac.Adec, and the multiplier that restores the
% fluid volume sum(phi) = Vcells is [eq:lambda]
%   lambda = (Vcells - sum(phi_a)) / sum(phi_b).
%
% Box constraint 0 <= phi <= 1:
%   ac.vol_project = false : clip after the update (the runs of Sections
%                            5.1-5.3 of the paper); the volume is exact
%                            whenever clipping is inactive.
%   ac.vol_project = true  : (default) bound-constrained projection: lambda
%                            solves sum(min(max(phi_a + lambda*phi_b,0),1)) =
%                            Vcells by bisection (monotone since phi_b >= 0),
%                            so clipping does not perturb the volume. The
%                            closed-form lambda is kept when nothing is clipped.
% freeze: logical mask of cells held at phiFrozen (design-frozen band).
% ac.cvec: Dirichlet boundary contribution of L_h (2D cases), if present.
% info (per inner step): vol_err = sum(phi) - Vcells, phi_min / phi_max of the
% unclipped update, n_clip = number of clipped entries, lambda.

info = struct('vol_err', zeros(ac.Kin,1), 'phi_min', zeros(ac.Kin,1), ...
              'phi_max', zeros(ac.Kin,1), 'n_clip', zeros(ac.Kin,1), 'lambda', zeros(ac.Kin,1));
project = ~isfield(ac,'vol_project') || ac.vol_project;
ghat = sens / (norm(sens(:)) + 1e-14);
for it = 1:ac.Kin
    Fw = phi.*(1-phi).*(0.5-phi);
    G  = phi.^2.*(1-phi).^2;
    forcing = -Fw - ac.beta*G.*ghat;
    forcing(freeze) = 0;
    Gm = G; Gm(freeze) = 0;

    rhs0 = phi + ac.tau*forcing;
    rhs1 = ac.tau*Gm;
    if isfield(ac,'cvec')
        pa = ac.Adec \ (rhs0(:) + ac.cvec);
    else
        pa = ac.Adec \ rhs0(:);
    end
    pb = ac.Adec \ rhs1(:);
    pa = reshape(pa, size(phi)); pb = reshape(pb, size(phi));
    pa(freeze) = phiFrozen(freeze); pb(freeze) = 0;

    sb = sum(pb(:));
    if abs(sb) > 1e-12
        lambda = (Vcells - sum(pa(:))) / sb;      % exact for the unclipped update
        phiraw = pa + lambda*pb;
        if project && (any(phiraw(:) < 0) || any(phiraw(:) > 1))
            lambda = project_vol(pa(:), pb(:), Vcells, lambda);
            phiraw = pa + lambda*pb;
        end
    else
        lambda = 0;      % no interface yet (binary phi, G = 0): no multiplier
        phiraw = pa;
    end
    phi = min(max(phiraw,0),1);
    phi(freeze) = phiFrozen(freeze);
    info.vol_err(it) = sum(phi(:)) - Vcells;
    info.phi_min(it) = min(phiraw(:)); info.phi_max(it) = max(phiraw(:));
    info.n_clip(it)  = nnz(phiraw < 0 | phiraw > 1);
    info.lambda(it)  = lambda;
end
end

% ------------------------------------------------------------------
function lambda = project_vol(pa, pb, Vcells, lam0)
% Bisection for f(lambda) = sum(min(max(pa + lambda*pb, 0), 1)) - Vcells = 0,
% f nondecreasing since pb >= 0; lam0 (closed form) seeds the bracket. If the
% volume is unreachable under the box constraint, lam0 is kept.
f = @(l) sum(min(max(pa + l*pb, 0), 1)) - Vcells;
tol = 1e-12*max(Vcells, 1);
f0 = f(lam0);
if abs(f0) <= tol, lambda = lam0; return, end
step = max(abs(lam0), 1);
if f0 < 0
    lo = lam0; hi = lam0 + step;
    while f(hi) < 0 && step < 1e30, step = 2*step; hi = lam0 + step; end
else
    hi = lam0; lo = lam0 - step;
    while f(lo) > 0 && step < 1e30, step = 2*step; lo = lam0 - step; end
end
if f(hi) < 0 || f(lo) > 0, lambda = lam0; return, end
for k = 1:200
    mid = 0.5*(lo + hi); fm = f(mid);
    if abs(fm) <= tol || (hi - lo) <= 1e-14*max(1, abs(mid)), lambda = mid; return, end
    if fm < 0, lo = mid; else, hi = mid; end
end
lambda = 0.5*(lo + hi);
end
