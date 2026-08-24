function model = vi_wnl_model(config, operators)
%VI_WNL_MODEL Connect a full cylinder discretization to the WNL engine.
%
% This adapter deliberately does not use the reduced matrices named B/B0 in
% the present linear code. Those matrices couple acceleration harmonics and
% are not the descriptor mass matrix.
%
% Required operator callbacks
% ---------------------------
% operators.B0(spec)
%   Complete descriptor matrix. Identity appears in momentum rows and in
%   the interface kinematic row; algebraic rows contain zero.
%
% operators.Lhat(spec,k)
%   kth fast-time Fourier coefficient of the complete linear operator,
%   including bulk, incompressibility, interface, wall, and pinning rows.
%
% operators.C(a,b,specA,specB,specOut)
%   Spatially projected symmetric quadratic action for one pair of temporal
%   coefficients. It must include convection, ALE mapping, kinematic,
%   evaluation-at-interface, and traction geometry terms.
%
% operators.D(a,b,c,specA,specB,specC,specOut)
%   Corresponding symmetric cubic action.
%
% Optional field-level alternatives
% ---------------------------------
% operators.CField(fieldA,fieldB,specOut)
% operators.DField(fieldA,fieldB,fieldC,specOut)
%   Return all output temporal coefficients in one pseudospectral action.
%   These override the coefficient-by-coefficient callbacks and are much
%   faster for a complete ALE residual.
%
% Optional
% --------
% operators.P(spec,nOut,nIn,direction)
%   Parameter-detuning block, including -omega_2*B0*d_tau when frequency is
%   varied.
%
% operators.conjugateAdjointRows(conjugatedField,sourceSpec)
%   Optional equation-row transformation for a physically conjugated
%   algebraic adjoint. It is needed when the -m operator equals a row-
%   transformed conjugate of the +m operator rather than its entrywise
%   conjugate in the stored row ordering.
%
% operators.solveForced(spec,forcing,opts)
%   Optional model-specific nonresonant forced solve. It may return either
%   a vector or a structure with fields vector and diagnostics. The cylinder
%   implementation exploits the low-rank interface coupling among temporal
%   harmonics and is checked in the complete unequilibrated equations.

requiredConfig = {'omega', 'N', 'ndof'};
requiredOperators = {'B0', 'Lhat', 'C', 'D'};
for j = 1:numel(requiredConfig)
    if ~isfield(config, requiredConfig{j})
        error('vi_wnl_model:MissingConfig', ...
            'config.%s is required.', requiredConfig{j});
    end
end
for j = 1:numel(requiredOperators)
    if ~isfield(operators, requiredOperators{j}) || ...
            ~isa(operators.(requiredOperators{j}), 'function_handle')
        error('vi_wnl_model:MissingOperator', ...
            'operators.%s must be a function handle.', ...
            requiredOperators{j});
    end
end

builder = struct();
builder.omega = config.omega;
builder.N = config.N;
builder.ndof = config.ndof;
builder.mass = operators.B0;
builder.linearFourier = operators.Lhat;
builder.quadraticLocal = operators.C;
builder.cubicLocal = operators.D;
if isfield(operators, 'nonlinearFourierShifts')
    builder.quadraticFourierShifts = ...
        operators.nonlinearFourierShifts;
    builder.cubicFourierShifts = ...
        operators.nonlinearFourierShifts;
end
if isfield(operators, 'P') && isa(operators.P, 'function_handle')
    builder.detuningBlock = operators.P;
end
if isfield(operators,'solveForced') && ...
        isa(operators.solveForced,'function_handle')
    builder.blockSolve = operators.solveForced;
end
model = wnl_fourier_model(builder);
if isfield(operators, 'CField') && ...
        isa(operators.CField, 'function_handle')
    model.quadratic = operators.CField;
end
if isfield(operators, 'DField') && ...
        isa(operators.DField, 'function_handle')
    model.cubic = operators.DField;
end
if isfield(operators,'conjugateAdjointRows') && ...
        isa(operators.conjugateAdjointRows,'function_handle')
    model.conjugateAdjointRows = operators.conjugateAdjointRows;
end
model.viConfig = config;
model.viOperators = operators;
end
