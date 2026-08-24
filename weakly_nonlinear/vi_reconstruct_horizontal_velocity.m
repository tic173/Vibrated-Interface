function [uR, uTheta] = vi_reconstruct_horizontal_velocity(r, k, m, wZ)
%VI_RECONSTRUCT_HORIZONTAL_VELOCITY Poloidal velocity from incompressibility.
%
% For w(r,z)=wHat(z)*J_m(k*r)*exp(i*m*theta),
%
%   u_H = wHat_z/k^2 * grad_H[J_m(k*r)*exp(i*m*theta)].
%
% Inputs:
%   r   [Nr,1] radial coordinates
%   k   radial wavenumber
%   m   integer azimuthal wavenumber
%   wZ  [Nz,Nt] vertical derivative coefficients
%
% Outputs:
%   uR, uTheta [Nr,Nz,Nt]

r = r(:);
if k <= 0
    error('vi_reconstruct_horizontal_velocity:BadWavenumber', ...
        'k must be positive.');
end

gamma = besselj(abs(m), k * r);
gammaR = 0.5 * k * ...
    (besselj(abs(m) - 1, k * r) - ...
     besselj(abs(m) + 1, k * r));
gammaTheta = complex(zeros(size(r)));
nonzero = abs(r) > sqrt(eps);
gammaTheta(nonzero) = 1i * m * gamma(nonzero) ./ r(nonzero);

atAxis = ~nonzero;
if any(atAxis)
    if abs(m) == 1
        gammaR(atAxis) = k / 2.0;
        gammaTheta(atAxis) = 1i * m * k / 2.0;
    else
        gammaR(atAxis) = 0.0;
        gammaTheta(atAxis) = 0.0;
    end
end

[nr, ~] = size(r);
[nz, nt] = size(wZ);
uR = complex(zeros(nr, nz, nt));
uTheta = complex(zeros(nr, nz, nt));
for n = 1:nt
    uR(:, :, n) = (gammaR / k^2) * transpose(wZ(:, n));
    uTheta(:, :, n) = ...
        (gammaTheta / k^2) * transpose(wZ(:, n));
end
end
