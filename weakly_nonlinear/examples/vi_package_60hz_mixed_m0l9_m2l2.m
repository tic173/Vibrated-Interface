function output = vi_package_60hz_mixed_m0l9_m2l2()
%VI_PACKAGE_60HZ_MIXED_M0L9_M2L2 Assemble the certified mixed-mode cache.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'weakly_nonlinear'));
auditFile = fullfile(root,'weakly_nonlinear','data', ...
    ['vi_wnl_complete_60Hz_ag0-4_m0l9-m2l2_', ...
     'chebyshev_zppb21_radial_audit.mat']);
if ~isfile(auditFile)
    error('vi_package_60hz_mixed_m0l9_m2l2:MissingAudit', ...
        'The strict radial audit is missing: %s',auditFile);
end
savedAudit = load(vi_wnl_resolve_data_file(auditFile),'report');
report = savedAudit.report;
if ~isfield(report,'reliable') || ~report.reliable
    error('vi_package_60hz_mixed_m0l9_m2l2:UnreliableAudit', ...
        'The mixed-mode audit has not passed every reliability gate.');
end

highFile = report.files.g21{end};
savedHigh = load(vi_wnl_resolve_data_file(highFile),'output');
output = savedHigh.output;
wnl = output.weaklyNonlinear;
modes = wnl.modes(:);
scales = wnl_mode_amplitude_scale(modes);
columnScaleSquared = abs(scales(:).').^2;

gPhysical = report.certifiedGPhysicalPeak;
hPhysical = report.certifiedHPhysicalPeak;
wnl.g = gPhysical.*columnScaleSquared;
wnl.gPhysicalPeak = gPhysical;
wnl.phaseSensitiveG = hPhysical.*columnScaleSquared;
wnl.phaseSensitiveGPhysicalPeak = hPhysical;
wnl.linearCoefficients = report.certifiedLambda(:);
wnl.mu = report.certifiedMu(:);
wnl.amplitudeScalesToPeakZetaOverH = scales(:);
wnl.coefficientComputed = true(2);
wnl.forcedSolvesValid = true(2);
wnl.validCubicScaling = true(2);
wnl.forcedSolvesExploratoryUsable = true(2);
wnl.requestedCoefficientValidity = true;
wnl.requestedExploratoryCoefficientAvailability = true;
wnl.fullCoefficientMatrixComputed = true;
wnl.numericalCoefficientValidity = true;
wnl.exploratoryCoefficientAvailability = true;
wnl.quantitativelyValid = logical(wnl.slowEnvelopeValid);
wnl.smallAmplitudeTransientAvailable = true;
wnl.smallAmplitudeTransientExploratoryAvailable = true;
wnl.coefficientCertification = report;
wnl.coefficientAssembly = report.coefficientAssembly;
output.weaklyNonlinear = wnl;
output.codeRelease = 'V71-radially-certified-mixed-coefficient-assembly';

outputFile = fullfile(root,'weakly_nonlinear','data',vi_wnl_output_filename(output.input));
output.input.run.outputFile = outputFile;
output.input.run.recoveredModeFile = outputFile;
output.input.run.postprocessSavedCoefficientsOnly = true;
output.input.run.reuseRecoveredModes = true;
output.input.run.modeRecoveryOnly = false;
output.coefficientCertification = struct( ...
    'auditFile',auditFile, ...
    'highGridCoefficientFile',highFile, ...
    'report',report);

save(outputFile,'output','-v7');
fprintf('Saved certified mixed-mode cache %s\n',outputFile);
fprintf('G in physical peak-zeta/h coordinates:\n');
disp(gPhysical);
fprintf('H in physical peak-zeta/h coordinates:\n');
disp(hPhysical);
end
