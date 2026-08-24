function roots = bessel_derivative_root(order, count)
%BESSEL_DERIVATIVE_ROOT First positive roots of J_order'(x).
%
% The root at x=0 is excluded. This function supplies the dependency used
% by the cylindrical linear solvers in this repository.

validateattributes(order, {'numeric'}, ...
    {'real', 'finite', 'nonnegative', 'integer', 'scalar'});
validateattributes(count, {'numeric'}, ...
    {'real', 'finite', 'positive', 'integer', 'scalar'});

derivative = @(x) 0.5 * ...
    (besselj(order - 1, x) - besselj(order + 1, x));
roots = zeros(1, count);
step = pi / 8.0;
xLeft = max(1.0e-8, sqrt(eps));
fLeft = derivative(xLeft);
found = 0;
maxX = (count + order / 2 + 4) * pi;

while found < count && xLeft < maxX
    xRight = xLeft + step;
    fRight = derivative(xRight);
    if isfinite(fLeft) && isfinite(fRight) && fLeft * fRight < 0
        root = fzero(derivative, [xLeft, xRight]);
        if found == 0 || abs(root - roots(found)) > 1.0e-7
            found = found + 1;
            roots(found) = root;
        end
        xLeft = root + 1.0e-6;
        fLeft = derivative(xLeft);
    else
        xLeft = xRight;
        fLeft = fRight;
    end
end

if found < count
    error('bessel_derivative_root:RootSearchFailed', ...
        'Located only %d of %d requested roots.', found, count);
end
end
