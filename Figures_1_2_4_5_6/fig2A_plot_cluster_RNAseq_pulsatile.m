clear all
load clustered_RNAseqDan.mat

%% plot gene expression clusters
n_cluster = max(cl_RNAseqDan_p.cluster);
%define color scheme
colors = [[.69:.15:1 ; 0:.25:.6; 0:.2:.45]'
    [0:.2:.45; 0:.25:.6; .69:.15:1]'];
b=rgb('DarkSlateGray')
c=rgb('DimGray')
d=rgb('Gray')
%adding colors
colors=[colors([1:3],:);b;c;d;c;b]
Cluster_times = [0:9 12];
z_puls = zscore(cl_RNAseqDan_p.FCs(:,ismember(cl_RNAseqDan_p.time,Cluster_times)),...
    [],2);
%%
fig=figure('Renderer', 'painters', 'Position', [1 400 950 500 ])
suptitle('mRNA clustering of p53 targets under p53 pulsatile (n=297)')
for i=1:n_cluster
    subplot(1,n_cluster,i)
    plot(0:10,z_puls(cl_RNAseqDan_p.cluster==i,:)','color', colors(i,:),'linewidth',0.5);pbaspect([1 1 1])
    hold on
    plot(0:10,mean(z_puls(cl_RNAseqDan_p.cluster==i,:))','k','linewidth',2)
    %text(4,1,sprintf('cluster %.d: n=%.d',i,sum(cl_RNAseqDan_p.cluster==i)),'FontSize',7)
    title(sprintf('cluster %.d: n=%.d',i,sum(cl_RNAseqDan_p.cluster==i)),'FontSize',9)
    xlim([0 10])
    ylim([-2.5 2.5])
    set(gca,'xtick',0:10,'xticklabel',[0:9,12])
    xlabel('time[h]')
    ylabel('zscore')

end
xlabel('time [h]')
set(fig,'Units','Inches');
set(gcf,'color','w');
set(fig,'PaperOrientation','landscape');
print(fig,'fig2_RNAseq_pulsatile_5','-dpdf','-r0')
