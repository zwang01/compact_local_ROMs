function stop = stop_rule(cfg, k, dphi, Jhist)
% STOP_RULE  Stopping rule of the paper (Section 5, common settings): stop
% after iteration k > 5 when ||phi_{k+1} - phi_k||_inf < cfg.stop_tol and, if
% cfg.stop_dJ_tol > 0, the objective has changed by less than cfg.stop_dJ_tol
% relative over the preceding 50 iterations (so not before k = 50).
% cfg.stop_tol = 0 switches the rule off (fixed number of iterations).
stop = false;
if cfg.stop_tol > 0 && k > 5 && dphi < cfg.stop_tol
    if cfg.stop_dJ_tol <= 0
        stop = true;
    elseif k >= 50
        stop = abs(Jhist(k+1) - Jhist(k-49))/max(abs(Jhist(k+1)), 1e-30) < cfg.stop_dJ_tol;
    end
end
end
