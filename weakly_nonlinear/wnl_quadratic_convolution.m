function coeffOut = wnl_quadratic_convolution(x, y, specOut, ...
    localAction, fourierShifts)
%WNL_QUADRATIC_CONVOLUTION Temporal convolution of a local bilinear action.
%
% The legacy callback
%   localAction(a,b,specA,specB,specOut)
% represents a time-independent quadratic map.  A periodically modulated
% nonlinear boundary condition can instead use
%   localAction(a,b,specA,specB,specOut,k,nuA,nuB),
% where k is the Fourier shift supplied in fourierShifts and
% nuA=nA+sA, nuB=nB+sB.  This is needed when vibration multiplies nonlinear
% interface geometry.

if nargin < 5 || isempty(fourierShifts)
    fourierShifts = 0;
end
fourierShifts = fourierShifts(:).';
localArgumentCount = nargin(localAction);

coeffOut = complex(zeros(specOut.ndof, numel(specOut.n)));
for ia = 1:numel(x.spec.n)
    nuA = x.spec.n(ia) + x.spec.s;
    for ib = 1:numel(y.spec.n)
        nuB = y.spec.n(ib) + y.spec.s;
        for k = fourierShifts
            nOutReal = nuA + nuB + k - specOut.s;
            nOut = round(nOutReal);
            if abs(nOutReal - nOut) > 1.0e-10
                error('wnl_quadratic_convolution:FrequencyMismatch', ...
                    ['Quadratic output with Fourier shift %g did not ', ...
                     'map to an integer harmonic.'], k);
            end
            io = find(specOut.n == nOut, 1);
            if isempty(io)
                continue;
            end
            if localArgumentCount < 0 || localArgumentCount >= 8
                value = localAction(x.coeff(:, ia), ...
                    y.coeff(:, ib), x.spec, y.spec, specOut, ...
                    k, nuA, nuB);
            elseif k == 0
                value = localAction(x.coeff(:, ia), ...
                    y.coeff(:, ib), x.spec, y.spec, specOut);
            else
                continue;
            end
            coeffOut(:, io) = coeffOut(:, io) + value;
        end
    end
end
end
