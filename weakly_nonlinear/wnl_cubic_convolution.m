function coeffOut = wnl_cubic_convolution(x, y, z, specOut, ...
    localAction, fourierShifts)
%WNL_CUBIC_CONVOLUTION Temporal convolution of a local trilinear action.
%
% The extended callback
%   localAction(a,b,c,specA,specB,specC,specOut,k,nuA,nuB,nuC)
% supports a periodic Fourier shift k in nonlinear interface geometry.

if nargin < 6 || isempty(fourierShifts)
    fourierShifts = 0;
end
fourierShifts = fourierShifts(:).';
localArgumentCount = nargin(localAction);

coeffOut = complex(zeros(specOut.ndof, numel(specOut.n)));
for ia = 1:numel(x.spec.n)
    nuA = x.spec.n(ia) + x.spec.s;
    for ib = 1:numel(y.spec.n)
        nuB = y.spec.n(ib) + y.spec.s;
        for ic = 1:numel(z.spec.n)
            nuC = z.spec.n(ic) + z.spec.s;
            for k = fourierShifts
                nOutReal = nuA + nuB + nuC + k - specOut.s;
                nOut = round(nOutReal);
                if abs(nOutReal - nOut) > 1.0e-10
                    error('wnl_cubic_convolution:FrequencyMismatch', ...
                        ['Cubic output with Fourier shift %g did not ', ...
                         'map to an integer harmonic.'], k);
                end
                io = find(specOut.n == nOut, 1);
                if isempty(io)
                    continue;
                end
                if localArgumentCount < 0 || localArgumentCount >= 11
                    value = localAction(x.coeff(:, ia), ...
                        y.coeff(:, ib), z.coeff(:, ic), ...
                        x.spec, y.spec, z.spec, specOut, ...
                        k, nuA, nuB, nuC);
                elseif k == 0
                    value = localAction(x.coeff(:, ia), ...
                        y.coeff(:, ib), z.coeff(:, ic), ...
                        x.spec, y.spec, z.spec, specOut);
                else
                    continue;
                end
                coeffOut(:, io) = coeffOut(:, io) + value;
            end
        end
    end
end
end
