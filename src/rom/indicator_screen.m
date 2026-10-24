function [eta_s, eta_a, cal] = indicator_screen(cal, gate, etaE_s, etaE_a, rs, ra, dJdU, Kref, Unorm)
% INDICATOR_SCREEN  The indicators that the acceptance test of Algorithm 2
% compares with tau_eta.
%   gate = 'cal'  (the method of the paper): calibrated indicators [eq:etacal],
%                 eta~ = eta^E / c, with c_s, c_a the medians of eta^E/eta over
%                 the last cal.mem direct evaluations of the exact indicators
%                 [eq:etaP]. The first cal.warm reduced iterations after the
%                 warm-up use the exact indicators and seed the ratios.
%   gate = 'eucl' : the uncalibrated Euclidean residuals [eq:indicator_flow]
%                 (ablation of Appendix C).
% etaE_s, etaE_a are the Euclidean indicators of the reduced step. cal holds the
% ratio histories cal.rs, cal.ra and the constants cal.warm, cal.mem.
eta_s = etaE_s; eta_a = etaE_a;
if ~strcmp(gate, 'cal') || isempty(Kref), return, end
if numel(cal.rs) < cal.warm
    [eta_s, eta_a] = indicator_exact(rs, ra, dJdU, Kref, Unorm);
    cal.rs(end+1) = etaE_s/eta_s; cal.ra(end+1) = etaE_a/eta_a;
else
    eta_s = etaE_s / median(cal.rs(max(1,end-cal.mem+1):end));
    eta_a = etaE_a / median(cal.ra(max(1,end-cal.mem+1):end));
end
end
