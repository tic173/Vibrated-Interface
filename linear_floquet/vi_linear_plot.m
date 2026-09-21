function handles = vi_linear_plot(result)
%VI_LINEAR_PLOT Plot modal traces and the reconstructed interface history.
cfg=result.config;
handles(1)=figure('Color','w','Name','Linear Floquet modes','Position',[100 100 900 360]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
plot(result.timeSeconds,real(result.modalDisplacement_m).','LineWidth',1.2);
xlabel('t (s)'); ylabel('Real modal displacement (m)'); grid on;
labels=arrayfun(@(md) char(md.input.label),result.modes,'UniformOutput',false);
legend(labels,'Interpreter','none','Location','best');
nexttile;
if isfield(result,'sweep')
    T=result.sweep.table; [~,order]=sort(T.kh);
    plot(T.kh(order),T.growth_per_s(order),'o-','MarkerSize',3,'LineWidth',1.1); hold on;
    ix=result.sweep.selectedIndex;
    plot(T.kh(ix),T.growth_per_s(ix),'rp','MarkerSize',11,'MarkerFaceColor','r');
    yline(0,'k:'); xlabel('Effective wavenumber k_{eff}h');
    ylabel('Growth rate (s^{-1})'); grid on;
    title(sprintf('Selected k_{eff}h = %.4g',T.kh(ix)));
else
    bar([result.modes.growthRatePerSecond]); yline(0,'k:');
    xlabel('Requested spatial mode'); ylabel('Growth rate (s^{-1})');
end
if isempty(result.interface_m), return; end
handles(2)=figure('Color','w','Name','Linear interface dynamics','Position',[100 100 900 650]);
if strcmp(cfg.geometry.type,'cartesian2d')
    imagesc(result.grid.X,result.timeSeconds, ...
        reshape(result.interface_m,numel(result.grid.X),[]).');
    axis xy; xlabel('x (m)'); ylabel('t (s)'); colorbar;
    title('Interface displacement (m)');
else
    frame=unique(round(linspace(1,numel(result.timeSeconds),min(4,numel(result.timeSeconds)))));
    columns=min(2,numel(frame));
    tiledlayout(ceil(numel(frame)/columns),columns,'TileSpacing','compact');
    limit=max(abs(result.interface_m(:,:,frame)),[],'all');
    if limit==0, limit=eps; end
    for it=frame
        nexttile;
        surf(result.grid.X,result.grid.Y,result.interface_m(:,:,it),'EdgeColor','none');
        xlabel('x (m)'); ylabel('y (m)'); zlabel('\zeta (m)');
        title(sprintf('t = %.3g s',result.timeSeconds(it))); view(35,30); axis tight;
        zlim([-limit limit]); caxis([-limit limit]); pbaspect([1 1 .65]);
    end
end
end
