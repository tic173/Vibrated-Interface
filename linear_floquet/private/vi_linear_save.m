function result = vi_linear_save(result)
base=char(result.config.output.directory);
if ~exist(base,'dir'), mkdir(base); end
[~,uniqueName]=fileparts(tempname(base));
folder=fullfile(base,['run_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '_' uniqueName]);
mkdir(folder); result.outputDirectory=folder;
mode=result.modes; input=[mode.input];
summary=table(string({input.label}).',[input.k].',[input.kStar].', ...
    [mode.growthRatePerSecond].',[mode.quasiFrequencyHz].', ...
    real([mode.multiplier].'),imag([mode.multiplier].'), ...
    'VariableNames',{'label','k_per_m','kh','growth_per_s','quasifrequency_Hz', ...
    'multiplier_real','multiplier_imag'});
if isfield(result,'sweep'), summary=result.sweep.table; end
writetable(summary,fullfile(folder,'growth_rates.csv'));
T=table(result.timeSeconds(:),'VariableNames',{'time_s'});
for i=1:numel(mode)
    T.(sprintf('mode%d_real_m',i))=real(result.modalDisplacement_m(i,:)).';
    T.(sprintf('mode%d_imag_m',i))=imag(result.modalDisplacement_m(i,:)).';
end
writetable(T,fullfile(folder,'modal_dynamics.csv'));
save(fullfile(folder,'linear_floquet.mat'),'result','-v7.3');
fprintf('Saved linear results to %s\n',folder);
end
