function cal = indicator_calibrate(cal, etaE_s, etaE_a, rs, ra, dJdU, Kref, Unorm, eta_s)
% INDICATOR_CALIBRATE  Update the calibration ratios c_s, c_a of [eq:etacal]
% at a full-order event of Algorithm 2 (a fallback or a verification), where
% the retained factorization is at hand: evaluate the exact indicators
% [eq:etaP] of the reduced pair just screened and append eta^E/eta. Also
% records cal.ratio = (indicator used by the gate)/(exact indicator), the
% effectivity of the calibrated state indicator at the event.
[ex_s, ex_a] = indicator_exact(rs, ra, dJdU, Kref, Unorm);
cal.rs(end+1) = etaE_s/ex_s; cal.ra(end+1) = etaE_a/ex_a;
cal.ratio(end+1) = eta_s/ex_s;
end
