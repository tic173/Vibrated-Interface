function [B0, L0, Lplus, Lminus, ...
    dLdaPlus, dLdaMinus, diagnostics] = ...
    vi_cylinder_wnl_linear_operators(discretization, m)
%VI_CYLINDER_WNL_LINEAR_OPERATORS Full primitive-variable cylinder block.
%
% The returned matrices satisfy
%   B0*q_t - L(tau)*q = 0,
%   L(tau)=L0+Lplus*exp(i*tau)+Lminus*exp(-i*tau).
% The flat-interface rows are the bulk momentum/divergence equations,
% rigid-wall and regularity rows, velocity and traction jumps, kinematic
% condition, and the selected free or pinned contact-line condition.

validateattributes(m, {'numeric'}, {'scalar', 'integer', 'finite'});
d = discretization;
p = d.parameters;
layout = d.layout;
nr = layout.nr;

[AD, BD, opsD] = layer_operator(nr, layout.nzD, d.r, ...
    d.Dr, d.Drr, d.DzD, d.DzzD, m, d.rhoD, d.muD*p.C);
[AL, BL, opsL] = layer_operator(nr, layout.nzL, d.r, ...
    d.Dr, d.Drr, d.DzL, d.DzzL, m, d.rhoL, d.muL*p.C);

ndof = layout.ndof;
A0 = sparse(ndof, ndof);
B0 = sparse(ndof, ndof);
A0([layout.d.ur, layout.d.ut, layout.d.w, layout.d.p], ...
    [layout.d.ur, layout.d.ut, layout.d.w, layout.d.p]) = AD;
B0([layout.d.ur, layout.d.ut, layout.d.w, layout.d.p], ...
    [layout.d.ur, layout.d.ut, layout.d.w, layout.d.p]) = BD;
A0([layout.l.ur, layout.l.ut, layout.l.w, layout.l.p], ...
    [layout.l.ur, layout.l.ut, layout.l.w, layout.l.p]) = AL;
B0([layout.l.ur, layout.l.ut, layout.l.w, layout.l.p], ...
    [layout.l.ur, layout.l.ut, layout.l.w, layout.l.p]) = BL;

Aplus = sparse(ndof, ndof);
Aminus = sparse(ndof, ndof);
AderivativePlus = sparse(ndof, ndof);
AderivativeMinus = sparse(ndof, ndof);

% No slip on the bottom and top, excluding the radial corner rows. The
% divergence row is retained at these horizontal boundaries.
for ir = 2:nr-1
    denseNode = grid_index(ir, 1, nr);
    lightNode = grid_index(ir, layout.nzL, nr);
    denseRows = [layout.d.ur(denseNode), ...
        layout.d.ut(denseNode), layout.d.w(denseNode)];
    lightRows = [layout.l.ur(lightNode), ...
        layout.l.ut(lightNode), layout.l.w(lightNode)];
    [A0, B0] = clear_rows(A0, B0, denseRows);
    [A0, B0] = clear_rows(A0, B0, lightRows);
    A0(denseRows(1), layout.d.ur(denseNode)) = 1;
    A0(denseRows(2), layout.d.ut(denseNode)) = 1;
    A0(denseRows(3), layout.d.w(denseNode)) = 1;
    A0(lightRows(1), layout.l.ur(lightNode)) = 1;
    A0(lightRows(2), layout.l.ut(lightNode)) = 1;
    A0(lightRows(3), layout.l.w(lightNode)) = 1;
end

% Axis regularity and the legacy free-slide sidewall conditions. These
% rows also cover the radial corners. For a Fourier mode m,
% u_r +/- i*u_theta have scalar orders |m +/- 1| at the axis.
[A0, B0] = impose_radial_boundaries(A0, B0, layout.d, ...
    layout.nzD, nr, d.r, d.Dr, m);
[A0, B0] = impose_radial_boundaries(A0, B0, layout.l, ...
    layout.nzL, nr, d.r, d.Dr, m);

% Join adjacent Chebyshev elements inside each fluid. Both copies of an
% artificial interface remain in the state vector. The left-copy rows
% impose continuity of velocity and pressure; the right-copy momentum rows
% impose continuity of the three vertical velocity derivatives. The
% right-copy pressure row remains the incompressibility equation. Because
% both elements use the same ALE map in one fluid, these smoothness
% conditions are linear even when the physical interface moves.
[A0, B0, denseMatchingRows] = impose_vertical_subdomain_matching( ...
    A0, B0, layout.d, nr, opsD.Dz, d.verticalD.interfacePairs);
[A0, B0, lightMatchingRows] = impose_vertical_subdomain_matching( ...
    A0, B0, layout.l, nr, opsL.Dz, d.verticalL.interfacePairs);

% Interface conditions are imposed at interior radial points. At the two
% radial corners, the regularity/sidewall rows above take precedence.
idD = layout.nzD;
idL = 1;
lapH = opsD.lapH;
for ir = 2:nr-1
    nodeD = grid_index(ir, idD, nr);
    nodeL = grid_index(ir, idL, nr);
    velocityRows = [layout.d.ur(nodeD), ...
        layout.d.ut(nodeD), layout.d.w(nodeD)];
    tractionRows = [layout.l.ur(nodeL), ...
        layout.l.ut(nodeL), layout.l.w(nodeL)];
    interfaceRows = [velocityRows, tractionRows];
    [A0, B0] = clear_rows(A0, B0, interfaceRows);
    Aplus(interfaceRows, :) = 0;
    Aminus(interfaceRows, :) = 0;
    AderivativePlus(interfaceRows, :) = 0;
    AderivativeMinus(interfaceRows, :) = 0;

    % Velocity continuity: lighter minus denser.
    A0(velocityRows(1), layout.l.ur(nodeL)) = 1;
    A0(velocityRows(1), layout.d.ur(nodeD)) = -1;
    A0(velocityRows(2), layout.l.ut(nodeL)) = 1;
    A0(velocityRows(2), layout.d.ut(nodeD)) = -1;
    A0(velocityRows(3), layout.l.w(nodeL)) = 1;
    A0(velocityRows(3), layout.d.w(nodeD)) = -1;

    rowD = nodeD;
    rowL = nodeL;
    radialDerivativeRow = d.Dr(ir, :);
    thetaDerivativeRow = sparse(1, ir, ...
        1i*m/max(d.r(ir), eps), 1, nr);

    % Tangential traction: T_d*n - T_l*n at the flat interface.
    A0(tractionRows(1), layout.d.ur) = ...
        d.muD*p.C*opsD.Dz(rowD, :);
    A0(tractionRows(1), layout.d.w( ...
        grid_index(1, idD, nr):grid_index(nr, idD, nr))) = ...
        A0(tractionRows(1), layout.d.w( ...
        grid_index(1, idD, nr):grid_index(nr, idD, nr))) + ...
        d.muD*p.C*radialDerivativeRow;
    A0(tractionRows(1), layout.l.ur) = ...
        -d.muL*p.C*opsL.Dz(rowL, :);
    A0(tractionRows(1), layout.l.w( ...
        grid_index(1, idL, nr):grid_index(nr, idL, nr))) = ...
        A0(tractionRows(1), layout.l.w( ...
        grid_index(1, idL, nr):grid_index(nr, idL, nr))) - ...
        d.muL*p.C*radialDerivativeRow;

    A0(tractionRows(2), layout.d.ut) = ...
        d.muD*p.C*opsD.Dz(rowD, :);
    A0(tractionRows(2), layout.d.w( ...
        grid_index(1, idD, nr):grid_index(nr, idD, nr))) = ...
        A0(tractionRows(2), layout.d.w( ...
        grid_index(1, idD, nr):grid_index(nr, idD, nr))) + ...
        d.muD*p.C*thetaDerivativeRow;
    A0(tractionRows(2), layout.l.ut) = ...
        -d.muL*p.C*opsL.Dz(rowL, :);
    A0(tractionRows(2), layout.l.w( ...
        grid_index(1, idL, nr):grid_index(nr, idL, nr))) = ...
        A0(tractionRows(2), layout.l.w( ...
        grid_index(1, idL, nr):grid_index(nr, idL, nr))) - ...
        d.muL*p.C*thetaDerivativeRow;

    % Normal traction, including static gravity and capillarity.
    normalRow = tractionRows(3);
    A0(normalRow, layout.l.p(nodeL)) = 1;
    A0(normalRow, layout.d.p(nodeD)) = -1;
    A0(normalRow, layout.l.w) = ...
        -2*d.muL*p.C*opsL.Dz(rowL, :);
    A0(normalRow, layout.d.w) = ...
        2*d.muD*p.C*opsD.Dz(rowD, :);
    A0(normalRow, layout.zeta) = ...
        d.hAtwood*p.g_sgn*sparse(1, ir, 1, 1, nr) - ...
        (1/p.Bd)*lapH(ir, :);

    gravityPlus = 0.5*p.aCritical*exp(1i*p.phase);
    gravityMinus = 0.5*p.aCritical*exp(-1i*p.phase);
    gravityDerivativePlus = 0.5*exp(1i*p.phase);
    gravityDerivativeMinus = 0.5*exp(-1i*p.phase);
    Aplus(normalRow, layout.zeta(ir)) = ...
        d.hAtwood*gravityPlus;
    Aminus(normalRow, layout.zeta(ir)) = ...
        d.hAtwood*gravityMinus;
    AderivativePlus(normalRow, layout.zeta(ir)) = ...
        d.hAtwood*gravityDerivativePlus;
    AderivativeMinus(normalRow, layout.zeta(ir)) = ...
        d.hAtwood*gravityDerivativeMinus;

    % Kinematic condition at the moving interface.
    kinematicRow = layout.zeta(ir);
    [A0, B0] = clear_rows(A0, B0, kinematicRow);
    A0(kinematicRow, layout.d.w(nodeD)) = -1;
    B0(kinematicRow, layout.zeta(ir)) = 1;
end

% Axis regularity and free/pinned contact-line condition for zeta.
zetaBoundaryRows = layout.zeta([1, nr]);
[A0, B0] = clear_rows(A0, B0, zetaBoundaryRows);
if m == 0
    A0(layout.zeta(1), layout.zeta) = d.Dr(1, :);
else
    A0(layout.zeta(1), layout.zeta(1)) = 1;
end
if strcmp(d.contactLine, 'pinned')
    A0(layout.zeta(nr), layout.zeta(nr)) = 1;
else
    A0(layout.zeta(nr), layout.zeta) = d.Dr(nr, :);
end

% Fixed volumes remove the spatially uniform m=0 interface displacement.
% The reduced linear code makes the same choice by retaining only positive
% Neumann Bessel roots. Replace one kinematic collocation row by the radial
% volume integral int_0^R zeta*r*dr=0.
volumeRow = [];
if m == 0
    volumeIr = volume_constraint_index(nr);
    volumeRow = layout.zeta(volumeIr);
    [A0, B0] = clear_rows(A0, B0, volumeRow);
    A0(volumeRow, layout.zeta) = ...
        trapezoidal_weights(d.r).'.*d.r.';
end

% Remove the common pressure gauge in the axisymmetric block.
gaugeRow = [];
if m == 0
    gaugeIr = max(2, min(nr-1, ceil(nr/2)));
    gaugeIz = d.denseGaugeVerticalIndex;
    gaugeNode = grid_index(gaugeIr, gaugeIz, nr);
    gaugeRow = layout.d.p(gaugeNode);
    [A0, B0] = clear_rows(A0, B0, gaugeRow);
    A0(gaugeRow, layout.d.p(gaugeNode)) = 1;
end

% wnl_fourier_model assembles i*omega*B0-Lhat. Since A0 is the
% zero-time-derivative residual, Lhat=-Ahat.
L0 = -A0;
Lplus = -Aplus;
Lminus = -Aminus;
dLdaPlus = -AderivativePlus;
dLdaMinus = -AderivativeMinus;

diagnostics = struct();
diagnostics.m = m;
diagnostics.gaugeRow = gaugeRow;
diagnostics.volumeRow = volumeRow;
diagnostics.contactLine = d.contactLine;
diagnostics.numberOfUnknowns = ndof;
diagnostics.massRankEstimate = sprank(B0);
diagnostics.sidewallConditions = d.sidewallConditions;
diagnostics.radialGridType = d.radial.typeUsed;
diagnostics.radialBasisConditionNumber = d.radial.conditionNumber;
diagnostics.denseVerticalElements = d.verticalD.numberOfElements;
diagnostics.lightVerticalElements = d.verticalL.numberOfElements;
diagnostics.denseMatchingRows = denseMatchingRows;
diagnostics.lightMatchingRows = lightMatchingRows;
diagnostics.numberOfVerticalMatchingRows = ...
    numel(denseMatchingRows)+numel(lightMatchingRows);
end

function [A, B, ops] = layer_operator(nr, nz, r, Dr, Drr, ...
    Dz, Dzz, m, rho, viscosity)
n = nr*nz;
ir = speye(nr);
iz = speye(nz);
rInverse = zeros(nr, 1);
rInverseSquared = zeros(nr, 1);
rInverse(2:end) = 1./r(2:end);
rInverseSquared(2:end) = 1./r(2:end).^2;
R1 = spdiags(rInverse, 0, nr, nr);
R2 = spdiags(rInverseSquared, 0, nr, nr);
lapH = Drr + R1*Dr - m^2*R2;
Dr2 = kron(iz, Dr);
Dz2 = kron(Dz, ir);
R12 = kron(iz, R1);
R22 = kron(iz, R2);
lap = kron(iz, lapH) + kron(Dzz, ir);

ur = 1:n;
ut = n+(1:n);
w = 2*n+(1:n);
pressure = 3*n+(1:n);
A = sparse(4*n, 4*n);
B = sparse(4*n, 4*n);

A(ur, ur) = -viscosity*(lap-R22);
A(ur, ut) = 2i*m*viscosity*R22;
A(ur, pressure) = Dr2;
A(ut, ut) = -viscosity*(lap-R22);
A(ut, ur) = -2i*m*viscosity*R22;
A(ut, pressure) = 1i*m*R12;
A(w, w) = -viscosity*lap;
A(w, pressure) = Dz2;
A(pressure, ur) = Dr2+R12;
A(pressure, ut) = 1i*m*R12;
A(pressure, w) = Dz2;
B(ur, ur) = rho*speye(n);
B(ut, ut) = rho*speye(n);
B(w, w) = rho*speye(n);

ops = struct();
ops.Dr = Dr2;
ops.Dz = Dz2;
ops.lapH = lapH;
ops.lap = lap;
end

function [A, B] = impose_radial_boundaries(A, B, layer, nz, ...
    nr, r, Dr, m)
orderPlus = abs(m+1);
orderMinus = abs(m-1);
for iz = 1:nz
    axisNode = grid_index(1, iz, nr);
    sideNode = grid_index(nr, iz, nr);
    axisRows = [layer.ur(axisNode), layer.ut(axisNode), ...
        layer.w(axisNode), layer.p(axisNode)];
    sideRows = [layer.ur(sideNode), layer.ut(sideNode), ...
        layer.w(sideNode)];
    [A, B] = clear_rows(A, B, [axisRows, sideRows]);

    radialNodes = grid_index(1, iz, nr):grid_index(nr, iz, nr);
    if orderPlus == 0
        A(axisRows(1), layer.ur(radialNodes)) = Dr(1, :);
        A(axisRows(1), layer.ut(radialNodes)) = 1i*Dr(1, :);
    else
        A(axisRows(1), layer.ur(axisNode)) = 1;
        A(axisRows(1), layer.ut(axisNode)) = 1i;
    end
    if orderMinus == 0
        A(axisRows(2), layer.ur(radialNodes)) = Dr(1, :);
        A(axisRows(2), layer.ut(radialNodes)) = -1i*Dr(1, :);
    else
        A(axisRows(2), layer.ur(axisNode)) = 1;
        A(axisRows(2), layer.ut(axisNode)) = -1i;
    end
    if m == 0
        A(axisRows(3), layer.w(radialNodes)) = Dr(1, :);
        A(axisRows(4), layer.p(radialNodes)) = Dr(1, :);
    else
        A(axisRows(3), layer.w(axisNode)) = 1;
        A(axisRows(4), layer.p(axisNode)) = 1;
    end

    A(sideRows(1), layer.ur(sideNode)) = 1;
    A(sideRows(2), layer.ut(sideNode)) = 1;
    A(sideRows(2), layer.ut(radialNodes)) = ...
        A(sideRows(2), layer.ut(radialNodes)) + r(end)*Dr(end, :);
    A(sideRows(3), layer.w(radialNodes)) = Dr(end, :);
end
end

function [A, B, matchingRows] = ...
    impose_vertical_subdomain_matching(A, B, layer, nr, Dz, interfacePairs)
matchingRows = zeros(0, 1);
for interfaceIndex = 1:size(interfacePairs, 1)
    leftIz = interfacePairs(interfaceIndex, 1);
    rightIz = interfacePairs(interfaceIndex, 2);
    for ir = 2:nr-1
        leftNode = grid_index(ir, leftIz, nr);
        rightNode = grid_index(ir, rightIz, nr);
        valueRows = [layer.ur(leftNode), layer.ut(leftNode), ...
            layer.w(leftNode), layer.p(leftNode)];
        derivativeRows = [layer.ur(rightNode), ...
            layer.ut(rightNode), layer.w(rightNode)];
        rows = [valueRows, derivativeRows];
        [A, B] = clear_rows(A, B, rows);

        % C0 continuity of primitive variables.
        A(valueRows(1), layer.ur([leftNode, rightNode])) = [-1, 1];
        A(valueRows(2), layer.ut([leftNode, rightNode])) = [-1, 1];
        A(valueRows(3), layer.w([leftNode, rightNode])) = [-1, 1];
        A(valueRows(4), layer.p([leftNode, rightNode])) = [-1, 1];

        % C1 continuity of velocity. Dz is block diagonal, so this row is
        % the right-element derivative minus the left-element derivative.
        derivativeJump = Dz(rightNode, :)-Dz(leftNode, :);
        A(derivativeRows(1), layer.ur) = derivativeJump;
        A(derivativeRows(2), layer.ut) = derivativeJump;
        A(derivativeRows(3), layer.w) = derivativeJump;
        matchingRows = [matchingRows; rows(:)]; %#ok<AGROW>
    end
end
end

function node = grid_index(ir, iz, nr)
node = ir + (iz-1)*nr;
end

function [A, B] = clear_rows(A, B, rows)
A(rows, :) = 0;
B(rows, :) = 0;
end

function index = volume_constraint_index(nr)
index = max(2, min(nr-1, ceil(nr/2)));
end

function weights = trapezoidal_weights(x)
x = x(:);
weights = zeros(size(x));
weights(1) = (x(2)-x(1))/2;
weights(end) = (x(end)-x(end-1))/2;
weights(2:end-1) = (x(3:end)-x(1:end-2))/2;
end
