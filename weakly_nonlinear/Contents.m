% Weakly nonlinear Floquet analysis.
%
% Model construction
%   wnl_fourier_model              - Assemble a Fourier descriptor model.
%   vi_wnl_model                   - Cylinder-operator adapter.
%   wnl_spec                       - Define one (m,s) Floquet block.
%   wnl_spec_lambda                - Read its continuous Floquet exponent.
%   cylinder_wnl_operators         - Full two-fluid cylinder WNL factory.
%   vi_cylinder_wnl_discretization - Grids, parameters, and state ordering.
%   vi_cylinder_wnl_linear_operators - Primitive linear descriptor blocks.
%
% Modes and forced fields
%   wnl_compute_mode               - Direct/left null vectors.
%   wnl_descriptor_residual        - Pressure-safe DAE mode residual.
%   wnl_complete_algebraic_state   - Minimum-norm pressure/constraint completion.
%   wnl_assert_mode_converged       - Gate coefficients on mode residuals.
%   wnl_reference_eigenvalue_consistency - Reduced/full exponent check.
%   wnl_assert_reference_eigenvalue_consistent - Gate model agreement.
%   wnl_conjugate_mode             - Physical conjugate Floquet mode.
%   wnl_conjugate_forced_solution  - Checked conjugate forced-field reuse.
%   wnl_solve_forced               - Bordered second-order solve.
%   vi_cylinder_wnl_forced_schur_solve - Low-rank cylinder temporal solve.
%   wnl_rank_aware_minimum_norm    - Equilibrated adaptive-rank forced solve.
%
% Coefficients
%   wnl_analyze_single_mode        - Detuning and self coefficient.
%   wnl_analyze_mode_set           - Self/cross coefficient matrix.
%   wnl_self_coefficient           - Mean/second-harmonic feedback.
%   wnl_cross_coefficient          - Cross-saturation coefficient.
%   wnl_find_quadratic_resonances  - Enumerate resonant triads.
%   wnl_rhs_landau                 - Coupled nonresonant amplitude ODE.
%
% Nonlinear actions
%   wnl_quadratic_convolution      - Two-field temporal convolution.
%   wnl_cubic_convolution          - Three-field temporal convolution.
%   wnl_fd_quadratic_action        - Residual-based directional C action.
%   wnl_fd_cubic_action            - Residual-based directional D action.
%   vi_cylinder_wnl_nonlinear_action - Cancellation-safe full ALE C/D actions.
%
% Vibrated-Interface helpers
%   vi_reconstruct_horizontal_velocity - Recover poloidal horizontal flow.
%   vi_dominant_floquet_root            - Multistart dominant-root search.
%   vi_operating_point_floquet_mode     - Exact mode at requested forcing.
%   vi_select_slow_operating_modes      - Select slow branches for two-mode WNL.
%   vi_compare_interface_dynamics       - Linear/WNL interface reconstruction.
%   vi_compare_interface_dynamics_modes - One/two-mode interface reconstruction.
%   vi_comparison_time_indices          - Select the transient comparison window.
%   vi_cubic_transient_correction       - Finite-time O(A^3) correction to linear growth.
%   vi_saved_small_amplitude_transient  - Fast saved-coefficient transient postprocessor.
%   vi_mode_initial_conditions          - Validate modal initial amplitudes/phases.
%   vi_resolve_analysis_amplitude       - Select an absolute forcing amplitude.
%   vi_cylinder_radial_grid             - Bessel-enriched radial operators.
%   vi_multidomain_chebyshev_grid       - Piecewise Chebyshev grid and maps.
%   vi_cylinder_vertical_grid           - Automatic/explicit cylinder z grid.
%   examples/vi_wnl_user_run_full_cylinder
%                                       - Physical WNL run and comparisons.
%   vi_wnl_apply_runtime_profile         - Development/balanced/final presets.
%   vi_wnl_apply_recovery_cache          - Validated direct/adjoint restart.
%   vi_save_output_record                - Robust absolute MAT-file saving.
%
% Verification
%   wnl_demo_stuart_landau         - Analytic coefficient test.
%   run_wnl_tests                  - Run all module tests.
