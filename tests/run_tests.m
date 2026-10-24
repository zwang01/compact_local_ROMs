% RUN_TESTS  All code-verification tests of Appendix C of the paper (a few
% minutes): finite-difference gradient checks, adjoint consistency, and the
% convergence orders of Table C.1 (Kovasznay, manufactured solution).
%   matlab -batch "setup_paths; run_tests"
t0 = tic;
fprintf('==================== test_gradient_fd ====================\n');        test_gradient_fd
fprintf('\n==================== test_adjoint_consistency ====================\n'); test_adjoint_consistency
fprintf('\n==================== test_kovasznay_order ====================\n');     test_kovasznay_order
fprintf('\n==================== test_mms_order ====================\n');           test_mms_order
fprintf('\nall tests finished in %.0f s\n', toc(t0));
