function [eta_s, eta_a] = indicator_exact(rs, ra, dJdU, Kref, Unorm)
% INDICATOR_EXACT  Preconditioned error indicators of the paper [eq:etaP]:
%   eta_s = ||K_ref^{-1} rs|| / ||U_ref||,
%   eta_a = ||K_ref^{-T} ra|| / ||K_ref^{-T} dJ/dU||,
% with K_ref the retained factorization (struct of solve closures, brk_lu.m)
% and Unorm the norm of the full-order state that goes with it. Three pairs of
% triangular solves.
eta_s = norm(Kref.solve(rs)) / Unorm;
eta_a = norm(Kref.solvet(ra)) / max(norm(Kref.solvet(dJdU)), 1e-30);
end
