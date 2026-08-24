function model = wnl_fourier_model(config)
%WNL_FOURIER_MODEL Build a WNL model from spatial Fourier coefficients.
%
% Required config fields
% ----------------------
% config.omega
%   Critical forcing frequency used in tau=omega*t.
%
% config.N
%   Temporal cutoff. Harmonic blocks use -N:N and subharmonic blocks use
%   -N-1:N.
%
% config.ndof
%   Number of spatial/constraint unknowns in one temporal harmonic.
%
% config.mass
%   Either the descriptor matrix B0 or a function B0=config.mass(spec).
%
% config.linearFourier
%   Lk=config.linearFourier(spec,k), where
%       L(tau)=sum_k Lk*exp(i*k*tau).
%   Return a sparse zero matrix for absent harmonics.
%
% Optional fields
% ---------------
% config.quadraticLocal
%   c=config.quadraticLocal(a,b,specA,specB,specOut).
%
% config.cubicLocal
%   d=config.cubicLocal(a,b,c,specA,specB,specC,specOut).
%
% config.detuningBlock
%   P=config.detuningBlock(spec,nOut,nIn,direction).
%   This callback may include -omega_2*B0*d_tau and other parameter
%   derivatives, so it is indexed by both input and output harmonics.
%
% config.blockSolve
%   Optional result=config.blockSolve(spec,forcing,opts) for a complete
%   nonresonant Floquet block. result may be a vector or a structure with
%   fields vector and diagnostics. Every caller still checks A*result-f.
%
% The assembled direct operator is
%   A_{n,n'} = [lambda+i*omega*(n+s)]*B0*delta_{n,n'} - L_{n-n'}.
% lambda=0 gives the original neutral-point formulation.  A nonzero lambda
% permits a shifted Floquet eigenproblem at a specified operating amplitude.
% The slow descriptor is kron(I,B0), without the factor omega.

required = {'omega', 'N', 'ndof', 'mass', 'linearFourier'};
for j = 1:numel(required)
    if ~isfield(config, required{j})
        error('wnl_fourier_model:MissingConfig', ...
            'config.%s is required.', required{j});
    end
end

model = struct();
model.config = config;
model.makeSpec = @(m, s, label) ...
    wnl_spec(m, s, config.N, config.ndof, label);
model.block = @(spec) assemble_block(config, spec);

if isfield(config, 'quadraticLocal') && ...
        isa(config.quadraticLocal, 'function_handle')
    quadraticShifts = 0;
    if isfield(config, 'quadraticFourierShifts')
        quadraticShifts = config.quadraticFourierShifts;
    end
    model.quadratic = @(x, y, specOut) ...
        wnl_quadratic_convolution(x, y, specOut, ...
        config.quadraticLocal, quadraticShifts);
end
if isfield(config, 'cubicLocal') && ...
        isa(config.cubicLocal, 'function_handle')
    cubicShifts = 0;
    if isfield(config, 'cubicFourierShifts')
        cubicShifts = config.cubicFourierShifts;
    end
    model.cubic = @(x, y, z, specOut) ...
        wnl_cubic_convolution(x, y, z, specOut, ...
        config.cubicLocal, cubicShifts);
end
if isfield(config, 'detuningBlock') && ...
        isa(config.detuningBlock, 'function_handle')
    model.detuning = @(field, direction, specOut) ...
        apply_detuning(config, field, direction, specOut);
end
end

function block = assemble_block(config, spec)
b0 = get_mass(config, spec);
if ~isequal(size(b0), [config.ndof, config.ndof])
    error('wnl_fourier_model:BadMass', ...
        'The descriptor matrix has the wrong size.');
end

nt = numel(spec.n);
lambda = wnl_spec_lambda(spec);
nTotal = config.ndof * nt;
A = sparse(nTotal, nTotal);
for row = 1:nt
    rowIdx = spatial_indices(row, config.ndof);
    for col = 1:nt
        colIdx = spatial_indices(col, config.ndof);
        k = spec.n(row) - spec.n(col);
        Lk = config.linearFourier(spec, k);
        if ~isequal(size(Lk), [config.ndof, config.ndof])
            error('wnl_fourier_model:BadLinearFourier', ...
                'linearFourier returned the wrong size for k=%d.', k);
        end
        value = -Lk;
        if row == col
            value = value + (lambda + ...
                1i * config.omega * (spec.n(row) + spec.s)) * b0;
        end
        if nnz(value) > 0
            A(rowIdx, colIdx) = value;
        end
    end
end

block = struct();
block.A = A;
block.Bslow = kron(speye(nt), sparse(b0));
block.spec = spec;
if isfield(config,'blockSolve') && isa(config.blockSolve,'function_handle')
    block.solve = @(forcing,requestedSpec,opts) ...
        config.blockSolve(requestedSpec,forcing,opts);
end
end

function coeffOut = apply_detuning(config, field, direction, specOut)
if ~wnl_equivalent_spec(field.spec, specOut)
    error('wnl_fourier_model:DetuningBlockMismatch', ...
        'Linear detuning must return to the input Floquet block.');
end

coeffOut = complex(zeros(specOut.ndof, numel(specOut.n)));
for row = 1:numel(specOut.n)
    for col = 1:numel(field.spec.n)
        pBlock = config.detuningBlock(specOut, ...
            specOut.n(row), field.spec.n(col), direction);
        if ~isequal(size(pBlock), [config.ndof, config.ndof])
            error('wnl_fourier_model:BadDetuningBlock', ...
                'detuningBlock returned the wrong matrix size.');
        end
        coeffOut(:, row) = coeffOut(:, row) + ...
            pBlock * field.coeff(:, col);
    end
end
end

function b0 = get_mass(config, spec)
if isa(config.mass, 'function_handle')
    b0 = config.mass(spec);
else
    b0 = config.mass;
end
end

function idx = spatial_indices(temporalIndex, ndof)
idx = (temporalIndex - 1) * ndof + (1:ndof);
end
