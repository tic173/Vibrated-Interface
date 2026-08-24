function [indices, actualEndPeriod] = vi_comparison_time_indices( ...
        forcingPeriods, requestedEndPeriod)
%VI_COMPARISON_TIME_INDICES Select a comparison record ending near a period.
%
% [indices,actualEndPeriod] = VI_COMPARISON_TIME_INDICES(periods,endPeriod)
% returns the contiguous indices from the start of the supplied time record
% through the sample nearest endPeriod. An empty endPeriod selects the whole
% record. The requested endpoint must lie inside the available record.

validateattributes(forcingPeriods, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonempty'});
forcingPeriods = forcingPeriods(:).';
if any(diff(forcingPeriods) < 0)
    error('vi_comparison_time_indices:NonmonotoneTime', ...
        'forcingPeriods must be nondecreasing.');
end

if nargin < 2 || isempty(requestedEndPeriod)
    endIndex = numel(forcingPeriods);
else
    validateattributes(requestedEndPeriod, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    timeTolerance = 64*eps(max(max(abs(forcingPeriods)), 1));
    if requestedEndPeriod < forcingPeriods(1)-timeTolerance || ...
            requestedEndPeriod > forcingPeriods(end)+timeTolerance
        error('vi_comparison_time_indices:EndOutsideRecord', ...
            ['Requested comparison endpoint %.12g is outside the ', ...
             'available forcing-period interval [%.12g, %.12g].'], ...
            requestedEndPeriod, forcingPeriods(1), forcingPeriods(end));
    end
    [~, endIndex] = min(abs(forcingPeriods-requestedEndPeriod));
end

indices = 1:endIndex;
actualEndPeriod = forcingPeriods(endIndex);
end
