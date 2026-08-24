function dxyz = wnl_fd_cubic_action(residual, q0, qTau0, ...
    x, xTau, y, yTau, z, zTau, tau, parameters, step)
%WNL_FD_CUBIC_ACTION Directional evaluation of D(x,y,z).
%
% For R(q)=A*q-C(q,q)-D(q,q,q)+..., D=-D^3R/6. The eight-point
% polarization below differentiates q and q_tau together.

if nargin < 13 || isempty(step)
    step = 2.0e-3;
end

d3 = zeros(size(residual(q0, qTau0, tau, parameters)), 'like', q0);
signValues = [-1, 1];
for sx = signValues
    for sy = signValues
        for sz = signValues
            q = q0 + step * (sx * x + sy * y + sz * z);
            qTau = qTau0 + step * ...
                (sx * xTau + sy * yTau + sz * zTau);
            d3 = d3 + sx * sy * sz * ...
                residual(q, qTau, tau, parameters);
        end
    end
end
d3 = d3 / (8.0 * step^3);
dxyz = -d3 / 6.0;
end
