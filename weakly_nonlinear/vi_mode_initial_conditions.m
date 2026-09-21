function initial = vi_mode_initial_conditions(modes, settings)
%VI_MODE_INITIAL_CONDITIONS Validate one- or two-mode Landau initial data.
%
% settings.amplitudesOverH contains nonnegative modal peak displacements and
% settings.phases contains their phases in radians. A self-conjugate
% axisymmetric H/SH mode has one real signed coordinate, so its phase must be
% 0 or pi modulo 2*pi. Distinct conjugate pairs retain a complex phase.

numberOfModes = numel(modes);
if numberOfModes < 1 || numberOfModes > 2
    error('vi_mode_initial_conditions:ModeCount', ...
        'The user driver supports either one or two retained Floquet modes.');
end
required = {'amplitudesOverH', 'phases'};
for fieldIndex = 1:numel(required)
    if ~isfield(settings, required{fieldIndex})
        error('vi_mode_initial_conditions:MissingInput', ...
            'input.initialConditions.%s is required.', ...
            required{fieldIndex});
    end
end

magnitudes = settings.amplitudesOverH(:);
phases = settings.phases(:);
validateattributes(magnitudes, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonnegative'});
validateattributes(phases, {'numeric'}, ...
    {'vector', 'real', 'finite'});
if numel(magnitudes) ~= numberOfModes || numel(phases) ~= numberOfModes
    error('vi_mode_initial_conditions:Size', ...
        ['Initial amplitudes and phases must each contain %d entries, ', ...
         'one for every retained mode.'], numberOfModes);
end
if magnitudes(1) <= 0
    error('vi_mode_initial_conditions:ZeroTarget', ...
        ['The first target mode must have a positive initial amplitude ', ...
         'for the linear/WNL comparison.']);
end

complexAmplitudes = magnitudes.*exp(1i*phases);
coordinateTypes = cell(numberOfModes,1);
labels = cell(numberOfModes, 1);
for modeIndex = 1:numberOfModes
    coordinateTypes{modeIndex} = ...
        wnl_mode_coordinate_type(modes(modeIndex));
    if strcmp(coordinateTypes{modeIndex},'real')
        phaseError = abs(sin(phases(modeIndex)));
        if phaseError > 1.0e-10
            error('vi_mode_initial_conditions:RealModePhase', ...
                ['Mode %d is self-conjugate and has a real signed ', ...
                 'amplitude. Its phase must be 0 or pi modulo 2*pi, ', ...
                 'not %.12g rad.'],modeIndex,phases(modeIndex));
        end
        signValue = sign(cos(phases(modeIndex)));
        if signValue == 0
            signValue = 1;
        end
        complexAmplitudes(modeIndex) = ...
            signValue*magnitudes(modeIndex);
    end
    if isfield(modes(modeIndex), 'label') && ...
            ~isempty(modes(modeIndex).label)
        labels{modeIndex} = char(modes(modeIndex).label);
    else
        labels{modeIndex} = sprintf('mode_%d', modeIndex);
    end
end

initial = struct();
initial.numberOfModes = numberOfModes;
initial.amplitudesOverH = magnitudes;
initial.phases = phases;
initial.complexAmplitudesOverH = complexAmplitudes;
initial.signedAmplitudesOverH = real(complexAmplitudes);
initial.coordinateTypes = coordinateTypes;
initial.labels = labels;
initial.zeroSeededModes = find(magnitudes == 0);
end
