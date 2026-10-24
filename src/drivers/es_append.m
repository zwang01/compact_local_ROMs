function es = es_append(es, k, etaE, etaEa, etaP, etaPa, etaC, etaCa, etaJ, err, erra, J, JF)
% ES_APPEND  One sample of the error study (cfg.errstudy = true): at an
% accepted reduced iteration k, the indicators available to Algorithm 2 and the
% true errors against a full-order solve at the same design.
%   etaE, etaEa : Euclidean state / adjoint residuals    [eq:indicator_flow]
%   etaP, etaPa : preconditioned (exact) indicators       [eq:etaP]
%   etaC, etaCa : the values the gate used (calibrated [eq:etacal] under 'cal')
%   etaJ        : dual-weighted estimate of J(U) - J(U_r) [eq:dwr_flow]
%   err, erra   : true relative state / adjoint errors
%   dJ          : true relative objective error |J(U_r) - J(U)|/|J(U)|
%   dJsigned    : J(U) - J(U_r), the quantity etaJ estimates
if isempty(es)
    es = struct('k',[],'etaE',[],'etaEa',[],'etaP',[],'etaPa',[],'etaC',[],'etaCa',[], ...
                'etaJ',[],'err',[],'erra',[],'dJ',[],'dJsigned',[]);
end
es.k(end+1) = k;
es.etaE(end+1) = etaE;   es.etaEa(end+1) = etaEa;
es.etaP(end+1) = etaP;   es.etaPa(end+1) = etaPa;
es.etaC(end+1) = etaC;   es.etaCa(end+1) = etaCa;
es.etaJ(end+1) = etaJ;
es.err(end+1) = err;     es.erra(end+1) = erra;
es.dJ(end+1) = abs(J - JF)/max(abs(JF), 1e-30);
es.dJsigned(end+1) = JF - J;
end
