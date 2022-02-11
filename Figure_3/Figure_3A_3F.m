clear all
load fit_results_all.mat
%% PLOT 1 GENE 
i=find(t_fit.gene=='BLOC1S2');
fig=figure('Renderer', 'painters', 'Position', [100 100 900 200])
%fig=figure('Renderer', 'painters', 'Position', [100 100 400 400])
%subplot(1,3,1)
%subplot(1,2,1)
subplot(1,4,1)
%subplot(2,2,1)
plot(0:9,t_fit.FC_RNA_p{i}(1:10), 'k','linewidth',2); hold on
title(sprintf('Input: %s mRNA dynamics', t_fit.gene(i)),'Fontsize',6)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
ylabel('FC')
legend(sprintf('%s exp-RNA Seq', t_fit.gene(i)),'Location', 'southeast')
legend boxoff 
%subplot(1,3,2)
%subplot(1,2,2)
subplot(1,4,2)
%subplot(2,2,2)
plot(0:9,t_fit.FC_Protein_p{i}(1:10),'ob','linewidth',1);
hold on
plot(0:9,t_fit.fitcurve_prot{i}, 'r','linewidth',2);
title(sprintf('Output: protein dynamics'),'Fontsize',6)
legend(sprintf('%s exp-Proteomics', t_fit.gene(i)),sprintf('%s fit', t_fit.gene(i)),'Location', 'southeast')
legend boxoff   
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
ylabel('FC')
txt=sprintf('R2 fit =%.2f',t_fit.r2_prot{i});
text(5,max(t_fit.fitcurve_prot{i})*0.9,txt,'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'baseline','Fontsize',6);
subplot(1,4,3)
%subplot(2,2,3)
plot(0:9,t_fit.FC_RNA_s{i}(1:10), 'k','linewidth',2); hold on
title(sprintf('Input: %s mRNA dynamics- p53 sustained', t_fit.gene(i)),'Fontsize',6)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
ylabel('FC')
legend(sprintf('%s exp-RNA Seq', t_fit.gene(i)),'Location', 'southeast')
legend boxoff 
subplot(1,4,4)
%subplot(2,2,4)
%subplot(1,3,3)
plot(0:9,t_fit.FC_Protein_s{i}(1:10),'ob','linewidth',1);
hold on
plot(0:9,t_fit.fitcurvepred_prot{i}, 'r','linewidth',2);
title(sprintf('Prediction: protein dynamics under sustained'),'Fontsize',6)
txt=sprintf('R2 pred =%.2f',t_fit.r2pred_prot{i});
text(5,max(t_fit.fitcurvepred_prot{i})*0.9,txt,'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'baseline','Fontsize',6);
legend('exp-Proteomics','prediction','Location', 'southeast')
legend boxoff 
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
ylabel('FC')
 
% subplot(1,3,3)
% plot(0:9,t_fit.FC_Protein_s{i}(1:10),'ob','linewidth',1);
% hold on
% plot(0:9,t_fit.fitcurvepred_prot{i}, 'r','linewidth',2);
% title(sprintf('Prediction: protein dynamics under sustained'))
% legend('exp-Proteomics','prediction','Location', 'southeast')
% legend boxoff 
% xticks([0:9])
% xticklabels({'0','1','2','3','4','5','6','7','8','9'})
% xlabel('Time (h)')
% ylabel('FC')
% txt=sprintf('R2=%.2f',t_fit.r2pred_prot{i});
% text(5,max(t_fit.fitcurvepred_prot{i}),txt)
set(gcf,'color','w');
orient(fig,'landscape')
%print('-bestfit','plot_general_MDM2_2A_longer','-dpdf')
print('-bestfit','plot_general_BLOC1S2_2A_longer','-dpdf')


%% PLOT 1 GENE FROM p53 to protein
i=find(t_fit.gene=='MDM2');
fig=figure('Renderer', 'painters', 'Position', [100 100 900 200])
% WB p53 values
p53_puls=[1 4.37 12.78 5.20 1.98 1.50 1.11 3.20 4.01 1.28 1.43 1.85 2.01];
p53_ps=smooth(p53_puls,3);
p53_sus=[1 4.37 12.78 15.74 25.69 22.41 21.50 36.88 25.02 37.43 29.85 15.41 22.63];
p53_ss=smooth(p53_sus,3);
subplot(1,6,1)
plot(0:9,p53_ps(1:10), 'k','linewidth',2); hold on
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xticks([0:9])
xlabel('Time (h)')
ylabel('FC')
title(sprintf('Input: p53 dynamics'),'Fontsize',6)
subplot(1,6,2)
plot(0:9,t_fit.FC_RNA_p{i}(1:10), 'ob','linewidth',1); hold on
plot(0:9,t_fit.fitcurve{i}(1:10),'k','linewidth',2);
title(sprintf('Input: %s mRNA dynamics', t_fit.gene(i)),'Fontsize',6)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
ylabel('FC')
legend(sprintf('%s exp-PNA Seq', t_fit.gene(i)),sprintf('%s fit', t_fit.gene(i)),'Location', 'southeast')
legend boxoff 
subplot(1,6,3)
plot(0:9,t_fit.FC_Protein_p{i}(1:10),'ob','linewidth',1);
hold on
% used to be 'r'
plot(0:9,t_fit.fitcurve_prot{i}, 'k','linewidth',2);
title(sprintf('Output: protein dynamics'),'Fontsize',6)
legend(sprintf('%s exp-Proteomics', t_fit.gene(i)),sprintf('%s fit', t_fit.gene(i)),'Location', 'southeast')
legend boxoff   
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
ylabel('FC')
txt=sprintf('R2 fit =%.2f',t_fit.r2_prot{i});
text(5,max(t_fit.fitcurve_prot{i})*0.9,txt,'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'baseline','Fontsize',6);
subplot(1,6,4)
% used to be 'r'
plot(0:9,t_fit.FC_RNA_s{i}(1:10), 'k','linewidth',2); hold on
title(sprintf('Input: %s mRNA dynamics- p53 sustained', t_fit.gene(i)),'Fontsize',6)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
ylabel('FC')
legend(sprintf('%s exp-RNA Seq', t_fit.gene(i)),'Location', 'southeast')
legend boxoff 
subplot(1,6,5)
plot(0:9,t_fit.FC_Protein_s{i}(1:10),'ob','linewidth',1);
hold on
% used to be 'r'
plot(0:9,t_fit.fitcurvepred_prot{i}, 'k','linewidth',2);
title(sprintf('Prediction: protein dynamics under sustained'),'Fontsize',6)
txt=sprintf('R2 pred =%.2f',t_fit.r2pred_prot{i});
text(5,max(t_fit.fitcurvepred_prot{i})*0.9,txt,'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'baseline','Fontsize',6);
legend('exp-Proteomics','prediction','Location', 'southeast')
legend boxoff 
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
ylabel('FC')
subplot(1,6,6)
plot(0:9,p53_ss(1:10), 'k','linewidth',2); hold on
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xticks([0:9])
xlabel('Time (h)')
ylabel('FC')
title(sprintf('p53 dynamics sustained'),'Fontsize',6)
set(gcf,'color','w');
orient(fig,'landscape')
%print('-bestfit','plot_general_MDM2_2A_longer','-dpdf')
print('-bestfit','plot_general_MDM2_2A_6pannels','-dpdf')


%%
% % plot figure type for the 17 genes, 4 at a time
% fig=figure('Renderer', 'painters', 'Position', [1 1 900 1500])
% set(gcf,'PaperType','A4') 
% %for i = 1 : 4
% %for i = 5 : 8
% %for i = 9 : 12
% %for i = 13 : 16
% %for i = 17 : 20
% for i = 21 : 24
% %for i = 25 : 28
% plotId=i-floor((i-1)/4)*4
% subplot(4, 3, 3*(plotId-1)+1) ;
% plot(0:12,t_fit.FC_RNA_p{i}(1:13), 'k','linewidth',2); hold on
% %plot(0:9,t_fit.FC_RNA_p{i}(1:10), 'k','linewidth',2); hold on
% title(sprintf('Input: %s mRNA under puls. p53', t_fit.gene(i)),'fontsize',8)
% hold on 
% plot([9 9],[0 max(t_fit.FC_RNA_p{i})],'--','color','k')
% %xticks([0:9])
% %xticklabels({'0','1','2','3','4','5','6','7','8','9'})
% xticks([0:12])
% xlim([0,12])
% xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12'})
% xlabel('Time (h)')
% txt=sprintf('mRNA cluster %d',t_fit.cluster(i));
% text(5,max(t_fit.FC_RNA_p{i}(1:10)),txt,'FontSize', 9);
% ylabel('FC')
% subplot(4, 3, 3*(plotId-1)+2) ;
% plot(0:9,t_fit.FC_Protein_p{i}(1:10),'ob','linewidth',1);
% hold on
% plot(0:9,t_fit.fitcurve_prot{i}, 'r','linewidth',2);
% title(sprintf('Output: %s Protein under puls. p53', t_fit.gene(i)),'fontsize',8)
% legend(sprintf('%s exp', t_fit.gene(i)),sprintf('%s fit', t_fit.gene(i)),'Location', 'southeast')
% legend boxoff   
% xticks([0:9])
% xticklabels({'0','1','2','3','4','5','6','7','8','9'})
% xlabel('Time (h)')
% %ylabel('FC')
% txt=sprintf('corr=%.2f;kd=%.2f;kp=%.2f; R2=%.2f',t_fit.FC_Protein_corr_p{i},t_fit.kd_p{i},t_fit.kp_p{i},t_fit.r2_prot{i});
% pos=max(t_fit.fitcurve_prot{i});
% text(0.5,pos,txt,'FontSize', 9);
% subplot(4, 3, 3*(plotId-1)+3) ;
% plot(0:9,t_fit.FC_Protein_s{i}(1:10),'ob','linewidth',1);
% hold on
% plot(0:9,t_fit.fitcurvepred_prot{i}, 'r','linewidth',2);
% title(sprintf('%s Protein under sust. p53', t_fit.gene(i)))
% legend('Data','Model','Location', 'southeast')
% %lgd.FontSize = 8;
% legend boxoff 
% xticks([0:9])
% xticklabels({'0','1','2','3','4','5','6','7','8','9'})
% xlabel('Time (h)')
% %ylabel('FC')
% txt=sprintf('corr=%.2f; R2=%.2f',t_fit.FC_Protein_corr_s{i},t_fit.r2pred_prot{i});
% text(0.5,max(t_fit.fitcurvepred_prot{i}),txt,'FontSize', 9);
% end
% set(fig,'Units','Inches');
% set(gcf,'color','w');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print('-fillpage','Badge6','-dpdf')


%% PLOT 1 GENE 
i=find(t_fit.gene=='SESN1');
fig=figure('Renderer', 'painters', 'Position', [100 100 900 200])
%fig=figure('Renderer', 'painters', 'Position', [100 100 400 400])
%subplot(1,3,1)
%subplot(1,2,1)
subplot(1,4,1)
%subplot(2,2,1)
plot(0:9,t_fit.FC_RNA_p{i}(1:10), 'k','linewidth',2); hold on
title(sprintf('Input: %s mRNA dynamics', t_fit.gene(i)),'Fontsize',6)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
ylabel('FC')
legend(sprintf('%s exp-RNA Seq', t_fit.gene(i)),'Location', 'southeast')
legend boxoff 
%subplot(1,3,2)
%subplot(1,2,2)
subplot(1,4,2)
%subplot(2,2,2)
plot(0:9,t_fit.FC_Protein_p{i}(1:10),'ob','linewidth',1);
hold on
plot(0:9,t_fit.fitcurve_prot{i}, 'r','linewidth',2);
title(sprintf('Output: protein dynamics'),'Fontsize',6)
legend(sprintf('%s exp-Proteomics', t_fit.gene(i)),sprintf('%s fit', t_fit.gene(i)),'Location', 'southeast')
legend boxoff   
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
ylabel('FC')
ylim([1 1.15])
txt=sprintf('R2 fit =%.2f',t_fit.r2_prot{i});
text(5,max(t_fit.fitcurve_prot{i})*0.9,txt,'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'baseline','Fontsize',6);
subplot(1,4,3)
%subplot(2,2,3)
plot(0:9,t_fit.FC_RNA_s{i}(1:10), 'k','linewidth',2); hold on
title(sprintf('Input: %s mRNA dynamics- p53 sustained', t_fit.gene(i)),'Fontsize',6)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
ylabel('FC')
legend(sprintf('%s exp-RNA Seq', t_fit.gene(i)),'Location', 'southeast')
legend boxoff 
subplot(1,4,4)
%subplot(2,2,4)
%subplot(1,3,3)
plot(0:9,t_fit.FC_Protein_s{i}(1:10),'ob','linewidth',1);
hold on
plot(0:9,t_fit.fitcurvepred_prot{i}, 'r','linewidth',2);
title(sprintf('Prediction: protein dynamics under sustained'),'Fontsize',6)
txt=sprintf('R2 pred =%.2f',t_fit.r2pred_prot{i});
text(5,max(t_fit.fitcurvepred_prot{i})*0.9,txt,'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'baseline','Fontsize',6);
legend('exp-Proteomics','prediction','Location', 'southeast')
legend boxoff 
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
ylabel('FC')
 
% subplot(1,3,3)
% plot(0:9,t_fit.FC_Protein_s{i}(1:10),'ob','linewidth',1);
% hold on
% plot(0:9,t_fit.fitcurvepred_prot{i}, 'r','linewidth',2);
% title(sprintf('Prediction: protein dynamics under sustained'))
% legend('exp-Proteomics','prediction','Location', 'southeast')
% legend boxoff 
% xticks([0:9])
% xticklabels({'0','1','2','3','4','5','6','7','8','9'})
% xlabel('Time (h)')
% ylabel('FC')
% txt=sprintf('R2=%.2f',t_fit.r2pred_prot{i});
% text(5,max(t_fit.fitcurvepred_prot{i}),txt)
set(gcf,'color','w');
orient(fig,'landscape')
%print('-bestfit','plot_general_MDM2_2A_longer','-dpdf')
print('-bestfit','plot_general_SESN1_2A_longer','-dpdf')