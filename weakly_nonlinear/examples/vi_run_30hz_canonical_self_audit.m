function outputs = vi_run_30hz_canonical_self_audit(stage,families,nr)
%VI_RUN_30HZ_CANONICAL_SELF_AUDIT Certify mode-local 30 Hz self terms.
%
% vi_run_30hz_canonical_self_audit('recovery',[],18)
% vi_run_30hz_canonical_self_audit('coefficients',[],18)
% vi_run_30hz_canonical_self_audit('coefficients')
% report = vi_run_30hz_canonical_self_audit('verify')
%
% A linear exponent and a cubic self interaction belong to one Floquet
% mode, not to the companion retained in a two-mode reduction. This audit
% computes the m=0 and m=2 diagonal terms once on fixed family operators.
% The three two-mode cases then reuse those canonical diagonal entries;
% cross interactions remain pair-specific and are not changed here.

% The source files provide already validated direct/adjoint vectors. The
% coefficient run always re-solves the nonlinear forced fields with the
% current solver and writes a separate, release-tagged output.

% families may contain 'm0' and/or 'm2'.

% See also VI_RUN_30HZ_COMMON_OPERATOR_CASES.

if nargin < 1 || isempty(stage)
    stage = 'verify';
end
stage = validatestring(stage,{'recovery','coefficients','verify'});
if nargin < 2 || isempty(families)
    families = {'m0','m2'};
elseif ischar(families) || isstring(families)
    families = cellstr(families);
end
if nargin < 3 || isempty(nr)
    nr = 16;
end
validateattributes(nr,{'numeric'},{'scalar','integer','>=',11});
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
definitions = family_definitions(root,nr);
for index = 1:numel(families)
    families{index} = validatestring(families{index}, ...
        fieldnames(definitions));
end

if strcmp(stage,'verify')
    outputs = build_report(root,definitions,nr);
    return;
end

outputs = cell(numel(families),1);
for index = 1:numel(families)
    definition = definitions.(families{index});
    if strcmp(stage,'recovery')
        requiredFile = definition.templateFile;
    else
        requiredFile = definition.sourceFile;
    end
    if ~isfile(requiredFile)
        error('vi_run_30hz_canonical_self_audit:MissingSource', ...
            'Validated recovery source is missing: %s', ...
            requiredFile);
    end
    switch stage
        case 'recovery'
            saved = load(vi_wnl_resolve_data_file(definition.templateFile),'output');
            preferred = saved.output.operatorMetadata.radialGrid. ...
                selectedBasisLabels;
            input = recovery_input( ...
                saved.output.input,definition,nr,preferred);
        case 'coefficients'
            saved = load(vi_wnl_resolve_data_file(definition.sourceFile),'output');
            input = coefficient_input(saved.output.input,definition);
    end
    outputs{index} = vi_wnl_run_with_input(input);
end
end

function definitions = family_definitions(root,nr)
m0Template = fullfile(root,'weakly_nonlinear','data', ...
    'vi_wnl_commonradial_ag0-3_fHz-30_modes-m0l6-m0l2_Nr16.mat');
m2Template = fullfile(root,'weakly_nonlinear','data', ...
    'vi_wnl_ag0-3_fHz-30_modes-m2l6-m2l2.mat');
m0Output = fullfile(root,'weakly_nonlinear','data',sprintf( ...
    'vi_wnl_canonicalself_ag0-3_fHz-30_modes-m0l6-m0l2_Nr%d.mat',nr));
m2Output = fullfile(root,'weakly_nonlinear','data',sprintf( ...
    'vi_wnl_canonicalself_ag0-3_fHz-30_modes-m2l6-m2l2_Nr%d.mat',nr));
definitions.m0 = struct( ...
    'templateFile',m0Template,'sourceFile',m0Output, ...
    'outputFile',m0Output);
definitions.m2 = struct( ...
    'templateFile',m2Template,'sourceFile',m2Output, ...
    'outputFile',m2Output);
if nr == 16
    definitions.m0.sourceFile = m0Template;
    definitions.m2.sourceFile = m2Template;
end
end

function input = recovery_input(input,definition,nr,preferredBasisLabels)
input.execution.profile = 'final';
input.numerics.Nr = nr;
input.numerics.radialGrid.preferredBasisLabels = preferredBasisLabels;
input.options.coefficientPairs = [];
input.run.outputFile = definition.outputFile;
input.run.recoveredModeFile = definition.outputFile;
input.run.modeRecoveryOnly = true;
input.run.reuseRecoveredModes = false;
input.run.requireRecoveredModes = false;
input.run.allowModeRecoveryDuringCoefficientRun = false;
input.run.postprocessSavedCoefficientsOnly = false;
input.run.saveResults = true;
input.run.saveProfile = 'compact';
input.run.plotLinearTransient = false;
input.run.compareLinearAndWNL = false;
end

function input = coefficient_input(input,definition)
input.execution.profile = 'final';
input.execution.startParallelPool = true;
input.execution.parallelWorkers = 2;
input.options.coefficientPairs = [1,1;2,2];
input.options.forcedSolveResidualTolerance = 1.0e-6;
input.options.forcedExploratoryResidualTolerance = [];
input.options.forcedModelSolverOnly = true;
input.options.stopOnUnconvergedMode = true;
input.options.stopOnUnconvergedForcedSolve = true;
input.options.forcedFailFast = true;
input.options.realCoefficientImaginaryTolerance = 5.0e-6;
input.weaklyNonlinear.reference = 'analysisAmplitude';
input.weaklyNonlinear.transientModel = 'cubicEnvelope';
input.weaklyNonlinear.autoSelectSlowModes = false;
input.weaklyNonlinear.requireSlowModesForComparison = false;
input.comparison.runGrowthRateSweep = false;
input.comparison.plotInterfaceDynamics = false;
input.comparison.plotRadialModeCorrections = false;
input.comparison.plotSecondModeInitialAmplitudeSweep = false;
input.run.outputFile = definition.outputFile;
input.run.recoveredModeFile = definition.sourceFile;
input.run.saveResults = true;
input.run.saveProfile = 'compact';
input.run.weaklyNonlinear = true;
input.run.linearPreview = true;
input.run.plotLinearTransient = false;
input.run.plotHarmonicNorms = false;
input.run.compareLinearAndWNL = false;
input.run.postprocessSavedCoefficientsOnly = false;
input.run.modeRecoveryOnly = false;
input.run.reuseRecoveredModes = true;
input.run.requireRecoveredModes = true;
input.run.allowModeRecoveryDuringCoefficientRun = false;
end

function report = build_report(root,definitions,nr)
familyNames = fieldnames(definitions);
canonical = struct([]);
for familyIndex = 1:numel(familyNames)
    definition = definitions.(familyNames{familyIndex});
    if ~isfile(definition.outputFile)
        error('vi_run_30hz_canonical_self_audit:MissingResult', ...
            ['Canonical result is missing: %s\nRun the coefficient ', ...
             'stage first.'],definition.outputFile);
    end
    saved = load(vi_wnl_resolve_data_file(definition.outputFile),'output');
    output = saved.output;
    result = output.weaklyNonlinear;
    for modeIndex = 1:numel(result.modes)
        entry.label = result.modes{modeIndex}.spec.label;
        entry.family = familyNames{familyIndex};
        entry.sourceFile = definition.outputFile;
        entry.codeRelease = output.codeRelease;
        entry.lambda = result.linearCoefficients(modeIndex);
        entry.gSelf = result.g(modeIndex,modeIndex);
        entry.gSelfPhysicalPeak = ...
            result.gPhysicalPeak(modeIndex,modeIndex);
        entry.directResidual = result.modes{modeIndex}.directResidual;
        entry.adjointResidual = result.modes{modeIndex}.leftResidual;
        entry.forcedValid = result.self{modeIndex}.forcedSolvesValid && ...
            result.self{modeIndex}.validCubicScaling;
        if isempty(canonical)
            canonical = entry;
        else
            canonical(end+1,1) = entry; %#ok<AGROW>
        end
    end
end

caseLabels = {'m0l6-m0l2','m0l6-m2l2','m2l6-m2l2'};
caseModes = {{'m0_l6_subharmonic','m0_l2_harmonic'}, ...
    {'m0_l6_subharmonic','m2_l2_harmonic'}, ...
    {'m2_l6_subharmonic','m2_l2_harmonic'}};
cases = repmat(struct(),numel(caseLabels),1);
for caseIndex = 1:numel(caseLabels)
    cases(caseIndex).label = caseLabels{caseIndex};
    cases(caseIndex).modeLabels = caseModes{caseIndex};
    cases(caseIndex).lambda = complex(zeros(2,1));
    cases(caseIndex).gDiagonal = complex(zeros(2,1));
    cases(caseIndex).gPhysicalPeakDiagonal = complex(zeros(2,1));
    for modeIndex = 1:2
        selected = find(strcmp({canonical.label}, ...
            caseModes{caseIndex}{modeIndex}));
        if numel(selected) ~= 1
            error('vi_run_30hz_canonical_self_audit:ModeLookup', ...
                'Expected one canonical entry for %s.', ...
                caseModes{caseIndex}{modeIndex});
        end
        cases(caseIndex).lambda(modeIndex) = canonical(selected).lambda;
        cases(caseIndex).gDiagonal(modeIndex) = ...
            canonical(selected).gSelf;
        cases(caseIndex).gPhysicalPeakDiagonal(modeIndex) = ...
            canonical(selected).gSelfPhysicalPeak;
    end
end

shared(1).label = 'm0_l6_subharmonic';
shared(1).caseIndices = [1,2];
shared(2).label = 'm2_l2_harmonic';
shared(2).caseIndices = [2,3];
for index = 1:numel(shared)
    occurrence = find(strcmp(caseModes{shared(index).caseIndices(1)}, ...
        shared(index).label));
    a = shared(index).caseIndices(1);
    b = shared(index).caseIndices(2);
    occurrenceB = find(strcmp(caseModes{b},shared(index).label));
    shared(index).linearDifference = abs( ...
        cases(a).lambda(occurrence)-cases(b).lambda(occurrenceB));
    shared(index).selfDifference = abs( ...
        cases(a).gDiagonal(occurrence)-cases(b).gDiagonal(occurrenceB));
    shared(index).passed = shared(index).linearDifference == 0 && ...
        shared(index).selfDifference == 0;
end
report.schemaVersion = 'V1-canonical-mode-local-self-audit';
report.frequencyHz = 30;
report.accelerationOverG = 3;
report.Nr = nr;
report.canonicalModes = canonical;
report.cases = cases;
report.sharedModes = shared;
report.pairInvariancePassed = all([shared.passed]);
report.localSolverGatesPassed = all([canonical.forcedValid]);
report.passed = report.pairInvariancePassed && ...
    report.localSolverGatesPassed;
report.certificationScope = [ ...
    'Pair invariance and local mode/forced-field residual gates only; ', ...
    'radial convergence requires a separate Nr comparison.'];
filename = fullfile(root,'weakly_nonlinear','data',sprintf( ...
    'vi_wnl_30Hz_canonical_self_audit_Nr%d.mat',nr));
save(filename,'report','-v7');
fprintf('\nCanonical 30 Hz linear/self-coefficient audit, Nr=%d\n',nr);
for index = 1:numel(canonical)
    item = canonical(index);
    fprintf(['  %-22s lambda=%+.12g, g_self=%+.12g, ', ...
        'g_peak=%+.12g, mode residuals %.3e/%.3e, forced valid=%d\n'], ...
        item.label,real(item.lambda),real(item.gSelf), ...
        real(item.gSelfPhysicalPeak),item.directResidual, ...
        item.adjointResidual,item.forcedValid);
end
for index = 1:numel(shared)
    fprintf(['  shared %-22s |Delta lambda|=%.3e, ', ...
        '|Delta g_self|=%.3e, passed=%d\n'],shared(index).label, ...
        shared(index).linearDifference,shared(index).selfDifference, ...
        shared(index).passed);
end
fprintf(['  audit passed=%d (pair invariance + local solver gates; ', ...
    'radial convergence not implied)\nSaved %s\n\n'], ...
    report.passed,filename);
end
