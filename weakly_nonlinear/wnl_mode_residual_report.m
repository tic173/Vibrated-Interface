function report = wnl_mode_residual_report(mode, layout, numberToPrint)
%WNL_MODE_RESIDUAL_REPORT Split direct/adjoint residuals by row family.
%
% The primitive-variable matrix reuses momentum rows for wall and interface
% conditions.  The labels therefore identify storage/row families and list
% the kinds of equations that may occupy them.  The Euclidean contributions
% use the same descriptor-aware physical scaling as wnl_compute_mode, so
% their squares add to the square of the reported total residual (up to
% roundoff). Algebraic pressure/gauge variables are excluded only from the
% normalization denominator; every equation row remains in the numerator.

if nargin < 3 || isempty(numberToPrint)
    numberToPrint = 6;
end

A = mode.block.A;
directResidualVector = A*mode.vector;
adjointResidualVector = A'*mode.left;
numberOfHarmonics = numel(mode.spec.n);
ndof = mode.spec.ndof;
offsets = (0:numberOfHarmonics-1)*ndof;

names = { ...
    'dense ur rows (radial momentum/wall/interface velocity)', ...
    'dense utheta rows (azimuthal momentum/wall/interface velocity)', ...
    'dense w rows (vertical momentum/wall/interface velocity)', ...
    'dense p rows (incompressibility/gauge)', ...
    'light ur rows (radial momentum/wall/tangential traction)', ...
    'light utheta rows (azimuthal momentum/wall/tangential traction)', ...
    'light w rows (vertical momentum/wall/normal traction)', ...
    'light p rows (incompressibility)', ...
    'zeta rows (kinematic/contact-line/volume)'};
localFamilies = {layout.d.ur, layout.d.ut, layout.d.w, layout.d.p, ...
    layout.l.ur, layout.l.ut, layout.l.w, layout.l.p, layout.zeta};

[~,directDetails] = wnl_descriptor_residual( ...
    A,mode.block.Bslow,mode.vector,'direct');
[~,adjointDetails] = wnl_descriptor_residual( ...
    A,mode.block.Bslow,mode.left,'adjoint');
operatorScale = directDetails.operatorScale;
directDenominator = operatorScale*max( ...
    directDetails.physicalVectorNorm,eps);
adjointDenominator = operatorScale*max( ...
    adjointDetails.physicalVectorNorm,eps);
directNorm = max(norm(directResidualVector), eps);
adjointNorm = max(norm(adjointResidualVector), eps);

families = repmat(struct(), numel(names), 1);
for familyIndex = 1:numel(names)
    local = localFamilies{familyIndex}(:);
    stacked = reshape(local+offsets, [], 1);
    directByHarmonic = zeros(numberOfHarmonics, 1);
    adjointByHarmonic = zeros(numberOfHarmonics, 1);
    for harmonicPosition = 1:numberOfHarmonics
        indices = local+offsets(harmonicPosition);
        directByHarmonic(harmonicPosition) = ...
            norm(directResidualVector(indices));
        adjointByHarmonic(harmonicPosition) = ...
            norm(adjointResidualVector(indices));
    end
    [~, directPeak] = max(directByHarmonic);
    [~, adjointPeak] = max(adjointByHarmonic);
    families(familyIndex).name = names{familyIndex};
    families(familyIndex).directScaled = ...
        norm(directResidualVector(stacked))/directDenominator;
    families(familyIndex).directFraction = ...
        norm(directResidualVector(stacked))/directNorm;
    families(familyIndex).directPeakFrequency = ...
        mode.spec.n(directPeak)+mode.spec.s;
    families(familyIndex).adjointScaled = ...
        norm(adjointResidualVector(stacked))/adjointDenominator;
    families(familyIndex).adjointFraction = ...
        norm(adjointResidualVector(stacked))/adjointNorm;
    families(familyIndex).adjointPeakFrequency = ...
        mode.spec.n(adjointPeak)+mode.spec.s;
end

report = struct();
report.operatorScale = operatorScale;
report.directTotal = norm(directResidualVector)/directDenominator;
report.adjointTotal = norm(adjointResidualVector)/adjointDenominator;
report.directResidualDetails = directDetails;
report.adjointResidualDetails = adjointDetails;
report.families = families;

[~, directOrder] = sort([families.directScaled], 'descend');
[~, adjointOrder] = sort([families.adjointScaled], 'descend');
numberToPrint = min(numberToPrint, numel(families));

fprintf('Largest direct residual row families\n');
for rank = 1:numberToPrint
    family = families(directOrder(rank));
    fprintf('  %.3e  %5.1f%%  n+s=%+.3g  %s\n', ...
        family.directScaled, 100*family.directFraction, ...
        family.directPeakFrequency, family.name);
end
fprintf('Largest adjoint residual variable families\n');
for rank = 1:numberToPrint
    family = families(adjointOrder(rank));
    fprintf('  %.3e  %5.1f%%  n+s=%+.3g  %s\n', ...
        family.adjointScaled, 100*family.adjointFraction, ...
        family.adjointPeakFrequency, family.name);
end
end
