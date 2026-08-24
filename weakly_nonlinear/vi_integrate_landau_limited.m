function result = vi_integrate_landau_limited(timeRequested, A0, lambda, g, maximumAmplitude)
%VI_INTEGRATE_LANDAU_LIMITED Integrate only inside the WNL amplitude range.
% Stops cleanly when max_j |A_j| reaches maximumAmplitude.  Returned times
% are a prefix of timeRequested, so linear and nonlinear records stay aligned.

timeRequested = timeRequested(:);
if numel(timeRequested) < 2
    error('vi_integrate_landau_limited:TimeVector', ...
        'At least two requested times are required.');
end
validateattributes(maximumAmplitude, {'numeric'}, ...
    {'scalar','real','positive','finite'});
A0 = A0(:);
rhs = @(t,A) wnl_rhs_landau(t, A, lambda, g);
events = @(t,A) amplitude_event(t, A, maximumAmplitude);
odeOptions = odeset('RelTol', 1.0e-8, 'AbsTol', 1.0e-11, ...
    'Events', events);
[tRaw, ARaw, tEvent] = ode45(rhs, ...
    [timeRequested(1), timeRequested(end)], A0, odeOptions);

lastTime = tRaw(end);
keep = timeRequested <= lastTime + ...
    64*eps(max(1.0, abs(lastTime)));
time = timeRequested(keep);
if isempty(time)
    time = timeRequested(1);
end
amplitudes = interp1(tRaw, ARaw, time, 'pchip');
if isvector(amplitudes) && numel(A0) > 1
    amplitudes = reshape(amplitudes, [], numel(A0));
end

result = struct();
result.time = time;
result.amplitudes = amplitudes;
result.completedRequestedWindow = lastTime >= timeRequested(end) - ...
    64*eps(max(1.0, abs(timeRequested(end))));
result.stoppedAtAmplitudeLimit = ~isempty(tEvent);
result.stopTime = lastTime;
result.maximumAmplitude = maximumAmplitude;
result.maximumReachedAmplitude = max(abs(ARaw(:)));
end

function [value, isterminal, direction] = amplitude_event(~, A, limit)
value = limit - max(abs(A));
isterminal = 1;
direction = -1;
end
