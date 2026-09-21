function result = vi_linear_most_unstable(cfg)
%VI_LINEAR_MOST_UNSTABLE Scan every mode; reconstruct ONLY the largest growth.
% This is a sampled maximum, not a continuous optimization or a full-spectrum
% guarantee. A failed mode aborts the scan rather than silently omitting it.
growth=vi_linear_growth(cfg);
inputs=[growth.modes.input]; modes=growth.modes;
[maximum,index]=max([modes.growthRatePerSecond]);
summary=table((1:numel(modes)).',string({inputs.label}).', ...
    [inputs.k].',[inputs.kStar].',[modes.growthRatePerSecond].', ...
    [modes.quasiFrequencyHz].',real([modes.multiplier].'),imag([modes.multiplier].'), ...
    'VariableNames',{'mode_index','label','k_per_m','kh','growth_per_s', ...
    'quasifrequency_Hz','multiplier_real','multiplier_imag'});
selected=growth;
selected.sweep=struct('table',summary,'modes',modes,'selectedIndex',index, ...
    'selection','Largest real exponent among successfully checked roots at the requested sample points.');
selected.modes=modes(index);
selected.config.modes=cfg.modes(index);
selected.maximumGrowthRatePerSecond=maximum;
if maximum>0, description='most unstable'; else, description='least damped (no positive growth found)'; end
fprintf('Selected %s: kh=%.8g, k=%.8g 1/m, growth=%+.8g 1/s (%s).\n', ...
    inputs(index).label,inputs(index).kStar,inputs(index).k,maximum,description);
result=vi_linear_dynamics(selected);
end
