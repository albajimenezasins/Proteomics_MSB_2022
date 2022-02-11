
clear all
load meanProteomics_non_scaled.mat
load meanRNAseqDan.mat

% qPCRs 1 and 2 where also done in the new irradiator, but with no time
% shift, however still one of the two experiments showed a second peak and
% the other did not. Time shifts are important as seen on qPCRS 3 and 4,
% which we are going to chose. qPCR 3 and 4 have a 30 min shift, they even
% sometimes peak before the RNAseq..but we rather chose that than not seen
% the second peak. We could peak one from 1-2 and another from 3-4 but its not 
% so clean cause there is no active 30 min shift
% WB data was in the new irradiator also but no 30 min delay
qPCR1_puls=readtable('./Raw_Data/qPCR_validation_1.xlsx','Sheet','puls');
qPCR1_sust=readtable('./Raw_Data/qPCR_validation_1.xlsx','Sheet','sust');
qPCR2_puls=readtable('./Raw_Data/qPCR_validation_2.xlsx','Sheet','puls');
qPCR2_sust=readtable('./Raw_Data/qPCR_validation_2.xlsx','Sheet','sust');
qPCR3_puls=readtable('./Raw_Data/qPCR_validation_3_30mshift_newirrad.xlsx','Sheet','puls');
qPCR4_puls=readtable('./Raw_Data/qPCR_validation_4_30mshift_newirrad.xlsx','Sheet','puls');
WB1_puls=readtable('./Raw_Data/Protein Validation for Alba Replicate 1.xlsx','Sheet','puls');
WB2_puls=readtable('./Raw_Data/Western Quantification DDB2 and KLF4 (second attempt).xlsx','Sheet','puls');
WB3_puls=readtable('./Raw_Data/Protein Validation for Alba Replicate 2.xlsx','Sheet','puls');
WB4_puls=readtable('./Raw_Data/p21 Protein Western Quantification Replicate 2_pulsatile.xlsx','Sheet','puls');
WB1_sust=readtable('./Raw_Data/Protein Validation for Alba Replicate 1.xlsx','Sheet','sust');
WB2_sust=readtable('./Raw_Data/W-Blot Sustained Rep 2.xlsx','Sheet','sust');

%% everything together_maxed_out
FigH=figure;
FigH=figure('Renderer', 'painters', 'Position', [1 400 950 500 ])
set(gcf,'color','w');
extraInputs = {'fontsize',8}; % name, value pairs
%***********************%
subplot(2,5,1)
yyaxis left; ax = gca; ax.FontSize = 8;
%ax.ColorOrder = [rgb('DarkRed'); rgb('DarkRed'); rgb('Gray'); rgb('Gray'); rgb('Gray'); rgb('Gray')];
r_fc3_qPCR=qPCR3_puls.p53'
r_fc4_qPCR=qPCR4_puls.p53'
r_qPCR=[r_fc3_qPCR(1:10);r_fc4_qPCR(1:10)];
err_qPCR=std(r_qPCR);
h1=plot(0:9,mean(r_qPCR),'linewidth',1.5);
ylim([0 2])
hold on
h2=errorbar(0:9,mean(r_qPCR),err_qPCR);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
ylabel('qPCR FC',extraInputs{:})
yyaxis right
r_fc1=RNAseqDan_p.FC1(find(RNAseqDan_p.gene=='TP53'),[1:10,12]);
r_fc2=RNAseqDan_p.FC2(find(RNAseqDan_p.gene=='TP53'),[1:10,12]);
r=[r_fc1(1:10);r_fc2(1:10)];
err=std(r);
plot(0:9,mean(r),'linewidth',1.5);
ylim([0 2])
ylabel('RNASeq FC',extraInputs{:})
hold on
errorbar(0:9,mean(r),err);
%txt = sprintf('MDM2 RNA cl%.f',cluster)
title('p53 RNA','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]',extraInputs{:})
xlim([0 9])
%***********************%
subplot(2,5,6)
yyaxis left; ax = gca; ax.FontSize = 8;
wb_fc1=WB1_puls.p53'
wb_fc3=WB3_puls.p53'
w=[wb_fc1(1:10);wb_fc3(1:10)];
err=std(w);
h1=plot(0:9,mean(w),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(w),err);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
ylabel('WB FC',extraInputs{:})
yyaxis right
p_fc=Proteomics_p.FC(find(Proteomics_p.gene=='TP53'),:);
p_fc1=Proteomics_p.FC1(find(Proteomics_p.gene=='TP53'),:);
p_fc2=Proteomics_p.FC2(find(Proteomics_p.gene=='TP53'),:);
p=[p_fc1(1:10);p_fc2(1:10)]
err=std(p)
h1=plot(0:9,p_fc(1:10),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(p),err);
ylabel('MassSpec FC',extraInputs{:})
set(h1,'color',rgb('Green'));
set(h2,'color',rgb('Green'));
set(gca,'ycolor',rgb('Green'))
title('p53 Protein','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%***********************%
subplot(2,5,2)
yyaxis left; ax = gca; ax.FontSize = 8;
r_fc3_qPCR=qPCR3_puls.Mdm2'
r_fc4_qPCR=qPCR4_puls.Mdm2'
r_qPCR=[r_fc3_qPCR(1:10);r_fc4_qPCR(1:10)];
err_qPCR=std(r_qPCR);
h1=plot(0:9,mean(r_qPCR),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(r_qPCR),err_qPCR);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
ylabel('qPCR FC',extraInputs{:})
yyaxis right
r_fc1=RNAseqDan_p.FC1(find(RNAseqDan_p.gene=='MDM2'),[1:10,12]);
r_fc2=RNAseqDan_p.FC2(find(RNAseqDan_p.gene=='MDM2'),[1:10,12]);
r=[r_fc1(1:10);r_fc2(1:10)];
err=std(r);
plot(0:9,mean(r),'linewidth',1.5);
ylabel('RNASeq FC',extraInputs{:})
hold on
errorbar(0:9,mean(r),err);
%txt = sprintf('MDM2 RNA cl%.f',cluster)
title('Mdm2 RNA','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]',extraInputs{:})
xlim([0 9])
%***********************%
subplot(2,5,7)
yyaxis left; ax = gca; ax.FontSize = 8;
wb_fc1=WB1_puls.Mdm2'
wb_fc3=WB3_puls.Mdm2'
w=[wb_fc1(1:10);wb_fc3(1:10)];
err=std(w);
h1=plot(0:9,mean(w),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(w),err);
ylabel('WB FC',extraInputs{:})
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
yyaxis right
p_fc=Proteomics_p.FC(find(Proteomics_p.gene=='MDM2'),:);
p_fc1=Proteomics_p.FC1(find(Proteomics_p.gene=='MDM2'),:);
p_fc2=Proteomics_p.FC2(find(Proteomics_p.gene=='MDM2'),:);
p=[p_fc1(1:10);p_fc2(1:10)]
err=std(p)
h1=plot(0:9,p_fc(1:10),'linewidth',1.5);
hold on
ylabel('MassSpec FC',extraInputs{:})
h2=errorbar(0:9,mean(p),err);
set(h1,'color',rgb('Green'));
set(h2,'color',rgb('Green'));
set(gca,'ycolor',rgb('Green'))
title('Mdm2 Protein','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%***********************%
subplot(2,5,3)
yyaxis left; ax = gca; ax.FontSize = 8;
r_fc3_qPCR=qPCR3_puls.p21'
r_fc4_qPCR=qPCR4_puls.p21'
r_qPCR=[r_fc3_qPCR(1:10);r_fc4_qPCR(1:10)];
err_qPCR=std(r_qPCR);
h1=plot(0:9,mean(r_qPCR),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(r_qPCR),err_qPCR);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
ylabel('qPCR FC',extraInputs{:})
yyaxis right
r_fc1=RNAseqDan_p.FC1(find(RNAseqDan_p.gene=='CDKN1A'),[1:10,12]);
r_fc2=RNAseqDan_p.FC2(find(RNAseqDan_p.gene=='CDKN1A'),[1:10,12]);
r=[r_fc1(1:10);r_fc2(1:10)];
err=std(r);
plot(0:9,mean(r),'linewidth',1.5);
ylabel('RNASeq FC',extraInputs{:})
hold on
errorbar(0:9,mean(r),err);
title('p21 RNA','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%***********************%
subplot(2,5,8)
yyaxis left; ax = gca; ax.FontSize = 8;
wb_fc1=WB1_puls.p21'
% wb_fc3=WB3_puls.p21'
wb_fc4=WB4_puls.p21'
w=[wb_fc1(1:10);wb_fc4(1:10)];
err=std(w);
h1=plot(0:9,mean(w),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(w),err);
%OR
%wb_fc=WB1_puls.p21'
%plot(0:9,wb_fc(1:10),'linewidth',1.5);
ylabel('WB FC',extraInputs{:})
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
yyaxis right
p_fc=Proteomics_p.FC(find(Proteomics_p.gene=='CDKN1A'),:);
p_fc1=Proteomics_p.FC1(find(Proteomics_p.gene=='CDKN1A'),:);
p_fc2=Proteomics_p.FC2(find(Proteomics_p.gene=='CDKN1A'),:);
p=[p_fc1(1:10);p_fc2(1:10)]
err=std(p)
h1=plot(0:9,p_fc(1:10),'linewidth',1.5);
hold on
ylabel('MassSpec FC',extraInputs{:})
h2=errorbar(0:9,mean(p),err);
set(h1,'color',rgb('Green'));
set(h2,'color',rgb('Green'));
set(gca,'ycolor',rgb('Green'))
 title('p21 Protein','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%***********************%
subplot(2,5,4)
yyaxis left; ax = gca; ax.FontSize = 8;
r_fc3_qPCR=qPCR3_puls.Wip1'
r_fc4_qPCR=qPCR4_puls.Wip1'
r_qPCR=[r_fc3_qPCR(1:10);r_fc4_qPCR(1:10)];
err_qPCR=std(r_qPCR);
h1=plot(0:9,mean(r_qPCR),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(r_qPCR),err_qPCR);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
ylabel('qPCR FC',extraInputs{:})
yyaxis right
r_fc1=RNAseqDan_p.FC1(find(RNAseqDan_p.gene=='PPM1D'),[1:10,12]);
r_fc2=RNAseqDan_p.FC2(find(RNAseqDan_p.gene=='PPM1D'),[1:10,12]);
r=[r_fc1(1:10);r_fc2(1:10)];
err=std(r);
plot(0:9,mean(r),'linewidth',1.5);
ylabel('RNASeq FC',extraInputs{:})
hold on
errorbar(0:9,mean(r),err);
title('Wip1 RNA','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%***********************%
subplot(2,5,9)
yyaxis left; ax = gca; ax.FontSize = 8;
wb_fc1=WB1_puls.Wip1'
wb_fc3=WB3_puls.Wip1'
w=[wb_fc1(1:10);wb_fc3(1:10)];
err=std(w);
h1=plot(0:9,mean(w),'linewidth',1.5);
hold on
ylabel('WB FC',extraInputs{:})
h2=errorbar(0:9,mean(w),err);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray'))
yyaxis right
p_fc=Proteomics_p.FC(find(Proteomics_p.gene=='PPM1D'),:);
p_fc1=Proteomics_p.FC1(find(Proteomics_p.gene=='PPM1D'),:);
p_fc2=Proteomics_p.FC2(find(Proteomics_p.gene=='PPM1D'),:);
p=[p_fc1(1:10);p_fc2(1:10)]
err=std(p)
h1=plot(0:9,p_fc(1:10),'linewidth',1.5);
hold on
ylabel('MassSpec FC',extraInputs{:})
h2=errorbar(0:9,mean(p),err);
set(h1,'color',rgb('Green'));
set(h2,'color',rgb('Green'));
set(gca,'ycolor',rgb('Green'))
title('Wip1 Protein','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
subplot(2,5,5)
yyaxis left; ax = gca; ax.FontSize = 8;
r_fc3_qPCR=qPCR3_puls.DDB2'
r_fc4_qPCR=qPCR4_puls.DDB2'
r_qPCR=[r_fc3_qPCR(1:10);r_fc4_qPCR(1:10)];
err_qPCR=std(r_qPCR);
h1=plot(0:9,mean(r_qPCR),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(r_qPCR),err_qPCR);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
ylabel('qPCR FC',extraInputs{:})
yyaxis right
r_fc1=RNAseqDan_p.FC1(find(RNAseqDan_p.gene=='DDB2'),[1:10,12]);
r_fc2=RNAseqDan_p.FC2(find(RNAseqDan_p.gene=='DDB2'),[1:10,12]);
r=[r_fc1(1:10);r_fc2(1:10)];
err=std(r);
plot(0:9,mean(r),'linewidth',1.5);
ylabel('RNASeq FC',extraInputs{:})
title('Ddb2 RNA','FontSize',9)
hold on
errorbar(0:9,mean(r),err);
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%***********************%
subplot(2,5,10)
yyaxis left; ax = gca; ax.FontSize = 8;
wb_fc2=WB2_puls.DDB2'
wb_fc3=WB3_puls.DDB2'
w=[wb_fc2(1:10);wb_fc3(1:10)];
err=std(w);
h1=plot(0:9,mean(w),'linewidth',1.5);
hold on
ylabel('WB FC',extraInputs{:})
h2=errorbar(0:9,mean(w),err);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
yyaxis right
p_fc=Proteomics_p.FC(find(Proteomics_p.gene=='DDB2'),:);
p_fc1=Proteomics_p.FC1(find(Proteomics_p.gene=='DDB2'),:);
p_fc2=Proteomics_p.FC2(find(Proteomics_p.gene=='DDB2'),:);
p=[p_fc1(1:10);p_fc2(1:10)]
err=std(p)
h1=plot(0:9,p_fc(1:10),'linewidth',1.5);
hold on
ylabel('MassSpec FC',extraInputs{:})
h2=errorbar(0:9,mean(p),err);
set(h1,'color',rgb('Green'));
set(h2,'color',rgb('Green'));
set(gca,'ycolor',rgb('Green'))
title('Ddb2 Protein','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%%
set(FigH,'Units','Inches');
set(gcf,'color','w');
set(FigH,'PaperOrientation','landscape');
print(FigH,'fig1_controls_pulsatile','-dpdf','-r0')

%% SUSTAINED
FigH=figure;
FigH=figure('Renderer', 'painters', 'Position', [1 400 950 500 ])
set(gcf,'color','w');
extraInputs = {'fontsize',8}; % name, value pairs
%***********************%
subplot(2,5,1)
yyaxis left; ax = gca; ax.FontSize = 8;
%ax.ColorOrder = [rgb('DarkRed'); rgb('DarkRed'); rgb('Gray'); rgb('Gray'); rgb('Gray'); rgb('Gray')];
r_fc1_qPCR=qPCR1_sust.p53'
r_fc2_qPCR=qPCR2_sust.p53'
r_qPCR=[r_fc1_qPCR(1:10);r_fc2_qPCR(1:10)];
err_qPCR=std(r_qPCR);
h1=plot(0:9,mean(r_qPCR),'linewidth',1.5);
ylim([0 2])
hold on
h2=errorbar(0:9,mean(r_qPCR),err_qPCR);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
ylabel('qPCR FC',extraInputs{:})
yyaxis right
r_fc1=RNAseqDan_s.FC1(find(RNAseqDan_s.gene=='TP53'),[1:10,12]);
r_fc2=RNAseqDan_s.FC2(find(RNAseqDan_s.gene=='TP53'),[1:10,12]);
r=[r_fc1(1:10);r_fc2(1:10)];
err=std(r);
h1=plot(0:9,mean(r),'linewidth',1.5);
ylim([0 2])
ylabel('RNASeq FC',extraInputs{:})
hold on
h2=errorbar(0:9,mean(r),err);
set(h1,'color',rgb('MediumBlue'));
set(h2,'color',rgb('MediumBlue'));
set(gca,'ycolor',rgb('MediumBlue')) 
%txt = sprintf('MDM2 RNA cl%.f',cluster)
title('p53 RNA','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]',extraInputs{:})
xlim([0 9])
%***********************%
subplot(2,5,6)
yyaxis left; ax = gca; ax.FontSize = 8;
wb_fc1=WB1_sust.p53'
wb_fc3=WB2_sust.p53'
w=[wb_fc1(1:10);wb_fc3(1:10)];
err=std(w);
h1=plot(0:9,mean(w),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(w),err);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
ylabel('WB FC',extraInputs{:})
yyaxis right
p_fc=Proteomics_s.FC(find(Proteomics_s.gene=='TP53'),:);
p_fc1=Proteomics_s.FC1(find(Proteomics_s.gene=='TP53'),:);
p_fc2=Proteomics_s.FC2(find(Proteomics_s.gene=='TP53'),:);
p=[p_fc1(1:10);p_fc2(1:10)]
err=std(p)
h1=plot(0:9,p_fc(1:10),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(p),err);
ylabel('MassSpec FC',extraInputs{:})
set(h1,'color',rgb('Green'));
set(h2,'color',rgb('Green'));
set(gca,'ycolor',rgb('Green'))
title('p53 Protein','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%***********************%
subplot(2,5,2)
yyaxis left; ax = gca; ax.FontSize = 8;
r_fc1_qPCR=qPCR1_sust.Mdm2'
r_fc2_qPCR=qPCR2_sust.Mdm2'
r_qPCR=[r_fc1_qPCR(1:10);r_fc2_qPCR(1:10)];
err_qPCR=std(r_qPCR);
h1=plot(0:9,mean(r_qPCR),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(r_qPCR),err_qPCR);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
ylabel('qPCR FC',extraInputs{:})
yyaxis right
r_fc1=RNAseqDan_s.FC1(find(RNAseqDan_s.gene=='MDM2'),[1:10,12]);
r_fc2=RNAseqDan_s.FC2(find(RNAseqDan_s.gene=='MDM2'),[1:10,12]);
r=[r_fc1(1:10);r_fc2(1:10)];
err=std(r);
h1=plot(0:9,mean(r),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(r),err);
ylabel('RNASeq FC',extraInputs{:})
set(h1,'color',rgb('MediumBlue'));
set(h2,'color',rgb('MediumBlue'));
set(gca,'ycolor',rgb('MediumBlue')) 
%txt = sprintf('MDM2 RNA cl%.f',cluster)
title('Mdm2 RNA','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]',extraInputs{:})
xlim([0 9])
%***********************%
subplot(2,5,7)
yyaxis left; ax = gca; ax.FontSize = 8;
wb_fc1=WB1_sust.Mdm2'
wb_fc3=WB2_sust.MDM2'
w=[wb_fc1(1:10);wb_fc3(1:10)];
err=std(w);
h1=plot(0:9,mean(w),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(w),err);
ylim ([0,16])
ylabel('WB FC',extraInputs{:})
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
yyaxis right
p_fc=Proteomics_s.FC(find(Proteomics_s.gene=='MDM2'),:);
p_fc1=Proteomics_s.FC1(find(Proteomics_s.gene=='MDM2'),:);
p_fc2=Proteomics_s.FC2(find(Proteomics_s.gene=='MDM2'),:);
p=[p_fc1(1:10);p_fc2(1:10)]
err=std(p)
h1=plot(0:9,p_fc(1:10),'linewidth',1.5);
hold on
ylabel('MassSpec FC',extraInputs{:})
h2=errorbar(0:9,mean(p),err);
set(h1,'color',rgb('Green'));
set(h2,'color',rgb('Green'));
set(gca,'ycolor',rgb('Green'))
title('Mdm2 Protein','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%***********************%
subplot(2,5,3)
yyaxis left; ax = gca; ax.FontSize = 8;
r_fc1_qPCR=qPCR1_sust.p21'
r_fc2_qPCR=qPCR2_sust.p21'
r_qPCR=[r_fc1_qPCR(1:10);r_fc2_qPCR(1:10)];
err_qPCR=std(r_qPCR);
h1=plot(0:9,mean(r_qPCR),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(r_qPCR),err_qPCR);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
ylabel('qPCR FC',extraInputs{:})
yyaxis right
r_fc1=RNAseqDan_s.FC1(find(RNAseqDan_s.gene=='CDKN1A'),[1:10,12]);
r_fc2=RNAseqDan_s.FC2(find(RNAseqDan_s.gene=='CDKN1A'),[1:10,12]);
r=[r_fc1(1:10);r_fc2(1:10)];
err=std(r);
h1=plot(0:9,mean(r),'linewidth',1.5);
ylabel('RNASeq FC',extraInputs{:})
hold on
h2=errorbar(0:9,mean(r),err);
set(h1,'color',rgb('MediumBlue'));
set(h2,'color',rgb('MediumBlue'));
set(gca,'ycolor',rgb('MediumBlue')) 
title('p21 RNA','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%***********************%
subplot(2,5,8)
yyaxis left; ax = gca; ax.FontSize = 8;
wb_fc1=WB1_sust.p21'
wb_fc3=WB2_sust.p21'
w=[wb_fc1(1:10);wb_fc3(1:10)];
err=std(w);
h1=plot(0:9,mean(w),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(w),err);
ylabel('WB FC',extraInputs{:})
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
yyaxis right
p_fc=Proteomics_s.FC(find(Proteomics_s.gene=='CDKN1A'),:);
p_fc1=Proteomics_s.FC1(find(Proteomics_s.gene=='CDKN1A'),:);
p_fc2=Proteomics_s.FC2(find(Proteomics_s.gene=='CDKN1A'),:);
p=[p_fc1(1:10);p_fc2(1:10)]
err=std(p)
h1=plot(0:9,p_fc(1:10),'linewidth',1.5);
hold on
ylabel('MassSpec FC',extraInputs{:})
h2=errorbar(0:9,mean(p),err);
set(h1,'color',rgb('Green'));
set(h2,'color',rgb('Green'));
set(gca,'ycolor',rgb('Green'))
 title('p21 Protein','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%***********************%
subplot(2,5,4)
yyaxis left; ax = gca; ax.FontSize = 8;
r_fc1_qPCR=qPCR1_sust.Wip1'
r_fc2_qPCR=qPCR2_sust.Wip1'
r_qPCR=[r_fc1_qPCR(1:10);r_fc2_qPCR(1:10)];
err_qPCR=std(r_qPCR);
h1=plot(0:9,mean(r_qPCR),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(r_qPCR),err_qPCR);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
ylabel('qPCR FC',extraInputs{:})
yyaxis right
r_fc1=RNAseqDan_s.FC1(find(RNAseqDan_s.gene=='PPM1D'),[1:10,12]);
r_fc2=RNAseqDan_s.FC2(find(RNAseqDan_s.gene=='PPM1D'),[1:10,12]);
r=[r_fc1(1:10);r_fc2(1:10)];
err=std(r);
h1=plot(0:9,mean(r),'linewidth',1.5);
ylabel('RNASeq FC',extraInputs{:})
hold on
h2=errorbar(0:9,mean(r),err);
set(h1,'color',rgb('MediumBlue'));
set(h2,'color',rgb('MediumBlue'));
set(gca,'ycolor',rgb('MediumBlue'))
title('Wip1 RNA','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%***********************%
subplot(2,5,9)
yyaxis left; ax = gca; ax.FontSize = 8;
wb_fc1=WB1_sust.Wip1'
wb_fc3=WB2_sust.WIP1'
w=[wb_fc1(1:10);wb_fc3(1:10)];
err=std(w);
h1=plot(0:9,mean(w),'linewidth',1.5);
hold on
ylabel('WB FC',extraInputs{:})
h2=errorbar(0:9,mean(w),err);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray'))
yyaxis right
p_fc=Proteomics_s.FC(find(Proteomics_s.gene=='PPM1D'),:);
p_fc1=Proteomics_s.FC1(find(Proteomics_s.gene=='PPM1D'),:);
p_fc2=Proteomics_s.FC2(find(Proteomics_s.gene=='PPM1D'),:);
p=[p_fc1(1:10);p_fc2(1:10)]
err=std(p)
h1=plot(0:9,p_fc(1:10),'linewidth',1.5);
hold on
ylabel('MassSpec FC',extraInputs{:})
h2=errorbar(0:9,mean(p),err);
set(h1,'color',rgb('Green'));
set(h2,'color',rgb('Green'));
set(gca,'ycolor',rgb('Green'))
title('Wip1 Protein','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])

subplot(2,5,5)
yyaxis left; ax = gca; ax.FontSize = 8;
r_fc1_qPCR=qPCR1_sust.DDB2'
r_fc2_qPCR=qPCR2_sust.DDB2'
r_qPCR=[r_fc1_qPCR(1:10);r_fc2_qPCR(1:10)];
err_qPCR=std(r_qPCR);
h1=plot(0:9,mean(r_qPCR),'linewidth',1.5);
hold on
h2=errorbar(0:9,mean(r_qPCR),err_qPCR);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
ylabel('qPCR FC',extraInputs{:})
yyaxis right
r_fc1=RNAseqDan_s.FC1(find(RNAseqDan_s.gene=='DDB2'),[1:10,12]);
r_fc2=RNAseqDan_s.FC2(find(RNAseqDan_s.gene=='DDB2'),[1:10,12]);
r=[r_fc1(1:10);r_fc2(1:10)];
err=std(r);
h1=plot(0:9,mean(r),'linewidth',1.5);
ylabel('RNASeq FC',extraInputs{:})
title('Ddb2 RNA','FontSize',9)
hold on
h2=errorbar(0:9,mean(r),err);
set(h1,'color',rgb('MediumBlue'));
set(h2,'color',rgb('MediumBlue'));
set(gca,'ycolor',rgb('MediumBlue'))
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%***********************%
subplot(2,5,10)
yyaxis left; ax = gca; ax.FontSize = 8;
wb_fc1=WB2_sust.DDB2'
wb_fc2=WB2_sust.DDDB2'
w=[wb_fc1(1:10);wb_fc2(1:10)];
err=std(w);
h1=plot(0:9,mean(w),'linewidth',1.5);
%h1=plot(0:9,wb_fc1(1:10),'linewidth',1.5);
hold on
ylabel('WB FC',extraInputs{:})
h2=errorbar(0:9,mean(w),err);
set(h1,'color',rgb('Gray'));
set(h2,'color',rgb('Gray'));
set(gca,'ycolor',rgb('Gray')) 
yyaxis right
p_fc=Proteomics_s.FC(find(Proteomics_s.gene=='DDB2'),:);
p_fc1=Proteomics_s.FC1(find(Proteomics_s.gene=='DDB2'),:);
p_fc2=Proteomics_s.FC2(find(Proteomics_s.gene=='DDB2'),:);
p=[p_fc1(1:10);p_fc2(1:10)]
err=std(p)
h1=plot(0:9,p_fc(1:10),'linewidth',1.5);
hold on
ylabel('MassSpec FC',extraInputs{:})
h2=errorbar(0:9,mean(p),err);
set(h1,'color',rgb('Green'));
set(h2,'color',rgb('Green'));
set(gca,'ycolor',rgb('Green'))
title('Ddb2 Protein','FontSize',9)
xticks([0:9])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('time [hours]')
xlim([0 9])
%%
set(FigH,'Units','Inches');
set(gcf,'color','w');
set(FigH,'PaperOrientation','landscape');
print(FigH,'fig1_controls_sustained','-dpdf','-r0')




























