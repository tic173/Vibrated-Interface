function result = vi_wnl_driver(operators, critical)
%VI_WNL_DRIVER Run the WNL engine for a completed cylinder discretization.
%
% Required critical fields:
%   omega, N, ndof, m, s
%
% Optional:
%   label, direct, left, detuning, options
%
% The operator structure is documented in vi_wnl_model.m.

required = {'omega', 'N', 'ndof', 'm', 's'};
for j = 1:numel(required)
    if ~isfield(critical, required{j})
        error('vi_wnl_driver:MissingCriticalField', ...
            'critical.%s is required.', required{j});
    end
end

config = struct();
config.omega = critical.omega;
config.N = critical.N;
config.ndof = critical.ndof;
model = vi_wnl_model(config, operators);

label = sprintf('m%d_mode', critical.m);
if isfield(critical, 'label')
    label = critical.label;
end
spec = model.makeSpec(critical.m, critical.s, label);
if isfield(critical, 'direct')
    spec.direct = critical.direct;
end
if isfield(critical, 'left')
    spec.left = critical.left;
end

opts = struct();
if isfield(critical, 'options')
    opts = critical.options;
end
if isfield(critical, 'detuning')
    opts.detuning = critical.detuning;
end
result = wnl_analyze_single_mode(model, spec, opts);
end
