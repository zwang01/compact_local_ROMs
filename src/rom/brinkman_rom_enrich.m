function rom = brinkman_rom_enrich(rom, S, Usnap, Wsnap)
% BRINKMAN_ROM_ENRICH  Window update of Algorithm 2 after a full-order solve
% (warm-up, fallback or verification): append the state snapshot to the state
% window W^s and the adjoint snapshot to the adjoint window W^a, drop the
% oldest snapshots beyond rom.w_max, recompute the two POD bases
%   rom.Phi  = Psi_s (state,   rom.r  modes),
%   rom.PhiA = Psi_a (adjoint, rom.rA modes),
% each truncated at the relative singular value rom.pod_tol and capped at
% rom.r_max, and the design-independent products used by the reduced solves,
%   K(phi) *Psi_s = rom.K0Phi   + dvec.*Psi_s,
%   K(phi)'*Psi_a = rom.K0tPhiA + dvec.*Psi_a,
% the Brinkman diagonal dvec being the only design dependence.
% Snapshots are normalized to unit length before they enter a window.
rom.SnapU = push(rom.SnapU, Usnap, rom.w_max);
rom.SnapW = push(rom.SnapW, Wsnap, rom.w_max);
[rom.Phi,  rom.r ] = pod(rom.SnapU, rom.r_max, rom.pod_tol);
[rom.PhiA, rom.rA] = pod(rom.SnapW, rom.r_max, rom.pod_tol);
rom.K0Phi   = S.K0   * rom.Phi;
rom.K0tPhiA = S.K0.' * rom.PhiA;
end

% ------------------------------------------------------------------
function Snap = push(Snap, new, w_max)
nj = norm(new);
if nj > 0, new = new/nj; end
Snap = [Snap, new];
if size(Snap,2) > w_max, Snap = Snap(:, end-w_max+1:end); end
end

function [Phi, r] = pod(Snap, r_max, pod_tol)
[U0, Sv, ~] = svd(Snap, 'econ');
s = diag(Sv);
r = min(r_max, max(1, sum(s >= pod_tol*s(1))));
Phi = U0(:,1:r);
end
