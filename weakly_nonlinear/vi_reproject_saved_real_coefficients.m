function output = vi_reproject_saved_real_coefficients( ...
        filename,imaginaryTolerance,saveResult)
%VI_REPROJECT_SAVED_REAL_COEFFICIENTS Revalidate completed real coefficients.
% Use this only when all saved forced fields passed their strict equation
% residual gates and the raw coefficients were withheld solely by the
% real-coordinate imaginary-leakage diagnostic.

if nargin < 2 || isempty(imaginaryTolerance)
    imaginaryTolerance = 5.0e-6;
end
if nargin < 3 || isempty(saveResult)
    saveResult = true;
end
validateattributes(filename,{'char','string'},{'scalartext'});
validateattributes(imaginaryTolerance,{'numeric'}, ...
    {'scalar','real','positive','finite'});
validateattributes(saveResult,{'logical','numeric'},{'scalar'});
filename = char(filename);
saved = load(filename,'output');
if ~isfield(saved,'output') || ...
        ~isfield(saved.output,'weaklyNonlinear') || ...
        isempty(saved.output.weaklyNonlinear)
    error('vi_reproject_saved_real_coefficients:MissingResult', ...
        'The MAT file does not contain output.weaklyNonlinear.');
end

output = saved.output;
oldRelease = output.codeRelease;
options = struct('realCoefficientImaginaryTolerance', ...
    imaginaryTolerance);
output.weaklyNonlinear = wnl_reproject_real_mode_set( ...
    output.weaklyNonlinear,options);
output.input.options.realCoefficientImaginaryTolerance = ...
    imaginaryTolerance;
output.codeRelease = 'V65-real-coefficient-reprojection';
output.weaklyNonlinear.coefficientRealityReprojection.sourceRelease = ...
    oldRelease;
output.weaklyNonlinear.coefficientRealityReprojection.timestamp = ...
    datetime('now','TimeZone','local');
if isfield(output,'comparison')
    output.comparison.available = false;
    output.comparison.quantitativelyValid = false;
    output.comparison.status = ...
        'coefficients reprojected; rerun saved-coefficient trajectory';
end

fprintf('Reprojected real signed coefficient matrix from saved raw values:\n');
disp(output.weaklyNonlinear.g);
fprintf(['  numerical coefficient validity = %d\n', ...
    '  slow-envelope validity         = %d\n'], ...
    output.weaklyNonlinear.numericalCoefficientValidity, ...
    output.weaklyNonlinear.slowEnvelopeValid);
if logical(saveResult)
    save(filename,'output','-v7.3');
    fprintf('Updated %s without recomputing fluid fields.\n',filename);
end
end
