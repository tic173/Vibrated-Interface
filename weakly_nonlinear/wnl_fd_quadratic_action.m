function cxy = wnl_fd_quadratic_action(residual, q0, qTau0, ...
    x, xTau, y, yTau, tau, parameters, step)
%WNL_FD_QUADRATIC_ACTION Directional evaluation of C(x,y).
%
% residual(q,q_tau,tau,parameters) must be the complete fixed-domain
% residual with the convention
%   R(q)=A*q-C(q,q)-D(q,q,q)+...
%
% The same perturbation is applied to q and q_tau, so ALE terms containing
% time derivatives are differentiated consistently.

if nargin < 11 || isempty(step)
    step = 1.0e-4;
end

rpp = residual(q0 + step * (x + y), ...
    qTau0 + step * (xTau + yTau), tau, parameters);
rpm = residual(q0 + step * (x - y), ...
    qTau0 + step * (xTau - yTau), tau, parameters);
rmp = residual(q0 + step * (-x + y), ...
    qTau0 + step * (-xTau + yTau), tau, parameters);
rmm = residual(q0 - step * (x + y), ...
    qTau0 - step * (xTau + yTau), tau, parameters);

d2 = (rpp - rpm - rmp + rmm) / (4.0 * step^2);
cxy = -0.5 * d2;
end
