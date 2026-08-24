function [analysisAmplitude, details] = ...
        vi_resolve_analysis_amplitude(forcing, linearSettings, ...
        criticalAmplitude)
%VI_RESOLVE_ANALYSIS_AMPLITUDE Select the requested physical forcing level.
%
% The preferred input is forcing.analysisAmplitude. For compatibility with
% V15 and earlier input files, a missing or empty value falls back to
% criticalAmplitude+linearSettings.accelerationOffsetFromCritical. If both
% are absent, the neutral critical amplitude is used.

validateattributes(criticalAmplitude, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'});

hasExplicitAmplitude = isfield(forcing, 'analysisAmplitude') && ...
    ~isempty(forcing.analysisAmplitude);
hasLegacyOffset = isfield(linearSettings, ...
    'accelerationOffsetFromCritical') && ...
    ~isempty(linearSettings.accelerationOffsetFromCritical);

if hasExplicitAmplitude
    analysisAmplitude = forcing.analysisAmplitude;
    source = 'input.forcing.analysisAmplitude';
    usedLegacyOffset = false;
elseif hasLegacyOffset
    validateattributes( ...
        linearSettings.accelerationOffsetFromCritical, {'numeric'}, ...
        {'scalar', 'real', 'finite'});
    analysisAmplitude = criticalAmplitude + ...
        linearSettings.accelerationOffsetFromCritical;
    source = ['legacy input.linear.', ...
        'accelerationOffsetFromCritical'];
    usedLegacyOffset = true;
else
    analysisAmplitude = criticalAmplitude;
    source = 'critical amplitude (no separate analysis input)';
    usedLegacyOffset = false;
end

validateattributes(analysisAmplitude, {'numeric'}, ...
    {'scalar', 'real', 'finite'});
if analysisAmplitude < 0
    error('vi_resolve_analysis_amplitude:NegativeAmplitude', ...
        ['The vibration magnitude a/g0 must be nonnegative. Change the ', ...
         'forcing phase instead of using a negative amplitude.']);
end

details = struct();
details.analysisAmplitude = analysisAmplitude;
details.criticalAmplitude = criticalAmplitude;
details.detuningFromCritical = ...
    analysisAmplitude-criticalAmplitude;
details.source = source;
details.usedLegacyOffset = usedLegacyOffset;
details.ignoredLegacyOffset = ...
    hasExplicitAmplitude && hasLegacyOffset;
end
