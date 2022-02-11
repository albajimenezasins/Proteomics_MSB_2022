clear all
load fit_results_all.mat
%% MDM2
time = 0:9;
j=find(t_fit.gene=='MDM2');
% kd MDM2 0.3531 fitted with r2 0.9143
%t_fit.kd_prot{j}
%t_fit.r2_prot{j}
kd_MDM2=[0.01;0.03;0.06;0.10;0.20;t_fit.kd_prot{j};0.5;1];
kp_MDM2=repmat(t_fit.kp_prot{j}, 1, 8);
kdel_MDM2=repmat(t_fit.kdel_prot{j}, 1, 8);

for i=1:length(kd_MDM2)
RNA_fit_puls = fit([-1 0:9]', [1 t_fit.FC_RNA_p{j}(1:10)]','linearinterp');
%[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_MDM2(i),kd_MDM2(i),kdel_MDM2(i)], RNA_fit_puls,  t_fit.FC_Protein_s{j}(1:10), time);
[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_MDM2(i),kd_MDM2(i),kdel_MDM2(i)], RNA_fit_puls, kp_MDM2(i)*ones(length(time),1)/kd_MDM2(i), time);
t_fit_Fig4_protein.FC_varykd{i}=t_fit_Fig4_protein.c_varykd{i}./repmat(t_fit_Fig4_protein.c_varykd{i}(1),1,length(t_fit_Fig4_protein.c_varykd{i}))
t_fit_Fig4_protein.log2FC_varykd{i}=log2(t_fit_Fig4_protein.FC_varykd{i});
end
fig=figure('Renderer', 'painters', 'Position', [100 100 1800 300])
for i=1:length(kd_MDM2)
subplot(1,9,i)
if kd_MDM2(i)==t_fit.kd_prot{j}
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'linewidth',1,'color','k')
xticks([0:9])
xlim([0,9])
ylim([-0.4,2.5])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
txt=sprintf('kd fitted=%.2f',kd_MDM2(i));
text(2,2.2,txt,'FontSize', 9);
else
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'color','b')
xticks([0:9])
xlim([0,9])
ylim([-0.4,2.5])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
if i==1
ylabel('MDM2 protein log2(FC)')
end
txt=sprintf('kd=%.2f',kd_MDM2(i));
text(2,2.2,txt,'FontSize', 9);
end
end
set(gcf,'color','w');
orient(fig,'landscape')
print('-bestfit','plot_vary_kd_MDM2','-dpdf')
%% AEN
time = 0:9;
j=find(t_fit.gene=='AEN');
% kd AEN 0.4981 fitted with r2 0.9436
%t_fit.kd_prot{j}
%t_fit.r2_prot{j}
kd_AEN=[0.05;0.10;0.25;t_fit.kd_prot{j};0.6;0.75;1];
%kd_AEN=[0.01;0.03;0.06;0.10;0.30;t_fit.kd_prot{j};0.8;1];
kp_AEN=repmat(t_fit.kp_prot{j}, 1, 8);
kdel_AEN=repmat(t_fit.kdel_prot{j}, 1, 8);

for i=1:length(kd_AEN)
RNA_fit_puls = fit([-1 0:9]', [1 t_fit.FC_RNA_p{j}(1:10)]','linearinterp');
%[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_AEN(i),kd_AEN(i),kdel_AEN(i)], RNA_fit_puls,  t_fit.FC_Protein_s{j}(1:10), time);
[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_AEN(i),kd_AEN(i),kdel_AEN(i)], RNA_fit_puls, kp_AEN(i)*ones(length(time),1)/kd_AEN(i), time);
t_fit_Fig4_protein.FC_varykd{i}=t_fit_Fig4_protein.c_varykd{i}./repmat(t_fit_Fig4_protein.c_varykd{i}(1),1,length(t_fit_Fig4_protein.c_varykd{i}))
t_fit_Fig4_protein.log2FC_varykd{i}=log2(t_fit_Fig4_protein.FC_varykd{i});
end
fig=figure
set(fig, 'Position',[100 100 1800 300])
for i=1:length(kd_AEN)
subplot(1,9,i)
if kd_AEN(i)==t_fit.kd_prot{j}
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'linewidth',1,'color','k')
xticks([0:9])
xlim([0,9])
ylim([-0.4,2])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
txt=sprintf('kd fitted=%.2f',kd_AEN(i));
text(2,1.7,txt,'FontSize', 9);
else
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'color','b')
xticks([0:9])
xlim([0,9])
ylim([-0.4,2])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
if i==1
ylabel('AEN protein log2(FC)')
end
txt=sprintf('kd=%.2f',kd_AEN(i));
text(2,1.7,txt,'FontSize', 9);
end
end
set(gcf,'color','w');
orient(fig,'landscape')
print('-bestfit','plot_vary_kd_AEN','-dpdf')
%% XPC
time = 0:9;
j=find(t_fit.gene=='XPC');
% kd XPC 0.439 fitted with r2 0.9476
% kd DCP1B 0.0555 fitted with r2 0.9816
t_fit.kd_prot{j}
t_fit.r2_prot{j}

kd_XPC=[0.005;0.01;0.02;t_fit.kd_prot{j};0.08;0.4;0.8];
%kd_XPC=[0.005;0.01;0.02;0.03;t_fit.kd_prot{j};0.6;0.8;1];
kp_XPC=repmat(t_fit.kp_prot{j}, 1, 8);
kdel_XPC=repmat(t_fit.kdel_prot{j}, 1, 8);

for i=1:length(kd_XPC)
RNA_fit_puls = fit([-1 0:9]', [1 t_fit.FC_RNA_p{j}(1:10)]','linearinterp');
%[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_XPC(i),kd_XPC(i),kdel_XPC(i)], RNA_fit_puls,  t_fit.FC_Protein_s{j}(1:10), time);
[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_XPC(i),kd_XPC(i),kdel_XPC(i)], RNA_fit_puls, kp_XPC(i)*ones(length(time),1)/kd_XPC(i), time);
t_fit_Fig4_protein.FC_varykd{i}=t_fit_Fig4_protein.c_varykd{i}./repmat(t_fit_Fig4_protein.c_varykd{i}(1),1,length(t_fit_Fig4_protein.c_varykd{i}))
t_fit_Fig4_protein.log2FC_varykd{i}=log2(t_fit_Fig4_protein.FC_varykd{i});
end
fig=figure
set(fig, 'Position',[100 100 1800 300])
for i=1:length(kd_XPC)
subplot(1,9,i)
if kd_XPC(i)==t_fit.kd_prot{j}
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'linewidth',1,'color','k')
xticks([0:9])
xlim([0,9])
ylim([-0.4,2])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
txt=sprintf('kd fitted=%.2f',kd_XPC(i));
text(2,1.7,txt,'FontSize', 9);
else
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'color','b')
xticks([0:9])
xlim([0,9])
ylim([-0.4,2])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
if i==1
ylabel('XPC protein log2(FC)')
end
txt=sprintf('kd=%.2f',kd_XPC(i));
text(2,1.7,txt,'FontSize', 9);
end
end
set(gcf,'color','w');
orient(fig,'landscape')
print('-bestfit','plot_vary_kd_XPC','-dpdf')

%% DCP1B
time = 0:9;
j=find(t_fit.gene=='DCP1B');
% kd DCP1B 0.04
kd_DCP1B=[0.005;0.01;0.02;0.03;t_fit.kd_prot{j};0.6;0.8;1];
kp_DCP1B=repmat(t_fit.kp_prot{j}, 1, 9);
kdel_DCP1B=repmat(t_fit.kdel_prot{j}, 1, 9);

for i=1:length(kd_DCP1B)
RNA_fit_puls = fit([-1 0:9]', [1 t_fit.FC_RNA_p{j}(1:10)]','linearinterp');
%[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_DCP1B(i),kd_DCP1B(i),kdel_DCP1B(i)], RNA_fit_puls,  t_fit.FC_Protein_s{j}(1:10), time);
[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_DCP1B(i),kd_DCP1B(i),kdel_DCP1B(i)], RNA_fit_puls, kp_DCP1B(i)*ones(length(time),1)/kd_DCP1B(i), time);
t_fit_Fig4_protein.FC_varykd{i}=t_fit_Fig4_protein.c_varykd{i}./repmat(t_fit_Fig4_protein.c_varykd{i}(1),1,length(t_fit_Fig4_protein.c_varykd{i}))
t_fit_Fig4_protein.log2FC_varykd{i}=log2(t_fit_Fig4_protein.FC_varykd{i});
end
fig=figure
set(fig, 'Position',[100 100 1800 300])
for i=1:length(kd_DCP1B)
subplot(1,9,i)
if kd_DCP1B(i)==t_fit.kd_prot{j}
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'linewidth',1,'color','k')
xticks([0:9])
xlim([0,9])
ylim([-0.4,2])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
txt=sprintf('kd fitted=%.2f',kd_DCP1B(i));
text(2,1.7,txt,'FontSize', 9);
else
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'color','b')
xticks([0:9])
xlim([0,9])
ylim([-0.4,2])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
if i==1
ylabel('DCP1B protein log2(FC)')
end
txt=sprintf('kd=%.2f',kd_DCP1B(i));
text(2,1.7,txt,'FontSize', 9);
end
end
set(gcf,'color','w');
orient(fig,'landscape')
print('-bestfit','plot_vary_kd_DCP1B','-dpdf')

%% FDXR
time = 0:9;
j=find(t_fit.gene=='FDXR');
% kd FDXR 0.04
kd_FDXR=[0.005;0.01;0.02;t_fit.kd_prot{j};0.08;0.4;0.8];
kp_FDXR=repmat(t_fit.kp_prot{j}, 1, 9);
kdel_FDXR=repmat(t_fit.kdel_prot{j}, 1, 9);

for i=1:length(kd_FDXR)
RNA_fit_puls = fit([-1 0:9]', [1 t_fit.FC_RNA_p{j}(1:10)]','linearinterp');
%[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_FDXR(i),kd_FDXR(i),kdel_FDXR(i)], RNA_fit_puls,  t_fit.FC_Protein_s{j}(1:10), time);
[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_FDXR(i),kd_FDXR(i),kdel_FDXR(i)], RNA_fit_puls, kp_FDXR(i)*ones(length(time),1)/kd_FDXR(i), time);
t_fit_Fig4_protein.FC_varykd{i}=t_fit_Fig4_protein.c_varykd{i}./repmat(t_fit_Fig4_protein.c_varykd{i}(1),1,length(t_fit_Fig4_protein.c_varykd{i}))
t_fit_Fig4_protein.log2FC_varykd{i}=log2(t_fit_Fig4_protein.FC_varykd{i});
end
fig=figure
set(fig, 'Position',[100 100 1800 300])
for i=1:length(kd_FDXR)
subplot(1,9,i)
if kd_FDXR(i)==t_fit.kd_prot{j}
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'linewidth',1,'color','k')
xticks([0:9])
xlim([0,9])
ylim([-0.4,2])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
txt=sprintf('kd fitted=%.2f',kd_FDXR(i));
text(2,1.7,txt,'FontSize', 9);
else
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'color','b')
xticks([0:9])
xlim([0,9])
ylim([-0.4,2])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
if i==1
ylabel('FDXR protein log2(FC)')
end
txt=sprintf('kd=%.2f',kd_FDXR(i));
text(2,1.7,txt,'FontSize', 9);
end
end
set(gcf,'color','w');
orient(fig,'landscape')
print('-bestfit','plot_vary_kd_FDXR','-dpdf')
%% SYTL1
time = 0:9;
j=find(t_fit.gene=='SYTL1');
% kd SYTL1 0.0135
kd_SYTL1=[0.002;0.005;0.008;t_fit.kd_prot{j};0.2;0.3;0.6;1];
kp_SYTL1=repmat(t_fit.kp_prot{j}, 1, 9);
kdel_SYTL1=repmat(t_fit.kdel_prot{j}, 1, 9);

for i=1:length(kd_SYTL1)
RNA_fit_puls = fit([-1 0:9]', [1 t_fit.FC_RNA_p{j}(1:10)]','linearinterp');
%[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_SYTL1(i),kd_SYTL1(i),kdel_SYTL1(i)], RNA_fit_puls,  t_fit.FC_Protein_s{j}(1:10), time);
[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_SYTL1(i),kd_SYTL1(i),kdel_SYTL1(i)], RNA_fit_puls, kp_SYTL1(i)*ones(length(time),1)/kd_SYTL1(i), time);
t_fit_Fig4_protein.FC_varykd{i}=t_fit_Fig4_protein.c_varykd{i}./repmat(t_fit_Fig4_protein.c_varykd{i}(1),1,length(t_fit_Fig4_protein.c_varykd{i}))
t_fit_Fig4_protein.log2FC_varykd{i}=log2(t_fit_Fig4_protein.FC_varykd{i});
end
fig=figure
set(fig, 'Position',[100 100 1800 300])
for i=1:length(kd_SYTL1)
subplot(1,9,i)
if kd_SYTL1(i)==t_fit.kd_prot{j}
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'linewidth',1,'color','k')
xticks([0:9])
xlim([0,9])
if i>4 ylim([0,1.2])
else ylim([0,0.1])
end
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
txt=sprintf('kd fitted=%.2f',kd_SYTL1(i));
text(2,1.7,txt,'FontSize', 9);
else
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'color','b')
xticks([0:9])
if i>4 ylim([0,1.2])
else ylim([0,0.1])
end
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
if i==1
ylabel('SYTL1 protein log2(FC)')
end
txt=sprintf('kd=%.2f',kd_SYTL1(i));
text(2,1.7,txt,'FontSize', 9);
end
end
set(gcf,'color','w');
orient(fig,'landscape')
print('-bestfit','plot_vary_kd_SYTL1','-dpdf')
%% DDB2
time = 0:9;
j=find(t_fit.gene=='DDB2');
% kd DDB2 0.0253
kd_DDB2=[0.002;0.005;0.008;t_fit.kd_prot{j};0.2;0.3;0.6;1];
kp_DDB2=repmat(t_fit.kp_prot{j}, 1, 9);
kdel_DDB2=repmat(t_fit.kdel_prot{j}, 1, 9);

for i=1:length(kd_DDB2)
RNA_fit_puls = fit([-1 0:9]', [1 t_fit.FC_RNA_p{j}(1:10)]','linearinterp');
%[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_DDB2(i),kd_DDB2(i),kdel_DDB2(i)], RNA_fit_puls,  t_fit.FC_Protein_s{j}(1:10), time);
[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_DDB2(i),kd_DDB2(i),kdel_DDB2(i)], RNA_fit_puls, kp_DDB2(i)*ones(length(time),1)/kd_DDB2(i), time);
t_fit_Fig4_protein.FC_varykd{i}=t_fit_Fig4_protein.c_varykd{i}./repmat(t_fit_Fig4_protein.c_varykd{i}(1),1,length(t_fit_Fig4_protein.c_varykd{i}))
t_fit_Fig4_protein.log2FC_varykd{i}=log2(t_fit_Fig4_protein.FC_varykd{i});
end
fig=figure
set(fig, 'Position',[100 100 1800 300])
for i=1:length(kd_DDB2)
subplot(1,9,i)
if kd_DDB2(i)==t_fit.kd_prot{j}
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'linewidth',1,'color','k')
xticks([0:9])
xlim([0,9])
% if i>4 ylim([0,1.2])
% else ylim([0,0.1])
% end
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
txt=sprintf('kd fitted=%.2f',kd_DDB2(i));
text(2,1.7,txt,'FontSize', 9);
else
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'color','b')
xticks([0:9])
% if i>4 ylim([0,1.2])
% else ylim([0,0.1])
% end
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
if i==1
ylabel('DDB2 protein log2(FC)')
end
txt=sprintf('kd=%.2f',kd_DDB2(i));
text(2,1.7,txt,'FontSize', 9);
end
end
set(gcf,'color','w');
orient(fig,'landscape')
print('-bestfit','plot_vary_kd_DDB2','-dpdf')
%% DGKA
time = 0:9;
j=find(t_fit.gene=='DGKA');
% kd DGKA 0.0253
kd_DGKA=[0.002;0.005;0.008;t_fit.kd_prot{j};0.2;0.3;0.6;1];
kp_DGKA=repmat(t_fit.kp_prot{j}, 1, 9);
kdel_DGKA=repmat(t_fit.kdel_prot{j}, 1, 9);

for i=1:length(kd_DGKA)
RNA_fit_puls = fit([-1 0:9]', [1 t_fit.FC_RNA_p{j}(1:10)]','linearinterp');
%[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_DGKA(i),kd_DGKA(i),kdel_DGKA(i)], RNA_fit_puls,  t_fit.FC_Protein_s{j}(1:10), time);
[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_DGKA(i),kd_DGKA(i),kdel_DGKA(i)], RNA_fit_puls, kp_DGKA(i)*ones(length(time),1)/kd_DGKA(i), time);
t_fit_Fig4_protein.FC_varykd{i}=t_fit_Fig4_protein.c_varykd{i}./repmat(t_fit_Fig4_protein.c_varykd{i}(1),1,length(t_fit_Fig4_protein.c_varykd{i}))
t_fit_Fig4_protein.log2FC_varykd{i}=log2(t_fit_Fig4_protein.FC_varykd{i});
end
fig=figure
set(fig, 'Position',[100 100 1800 300])
for i=1:length(kd_DGKA)
subplot(1,9,i)
if kd_DGKA(i)==t_fit.kd_prot{j}
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'linewidth',1,'color','k')
xticks([0:9])
xlim([0,9])
% if i>4 ylim([0,1.2])
% else ylim([0,0.1])
% end
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
txt=sprintf('kd fitted=%.2f',kd_DGKA(i));
text(2,1.7,txt,'FontSize', 9);
else
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'color','b')
xticks([0:9])
% if i>4 ylim([0,1.2])
% else ylim([0,0.1])
% end
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
if i==1
ylabel('DGKA protein log2(FC)')
end
txt=sprintf('kd=%.2f',kd_DGKA(i));
text(2,1.7,txt,'FontSize', 9);
end
end
set(gcf,'color','w');
orient(fig,'landscape')
print('-bestfit','plot_vary_kd_DGKA','-dpdf')

%% TNFRSF10C 
time = 0:9;
j=find(t_fit.gene=='TNFRSF10C');
% kd TNFRSF10C 0.3365
kd_TNFRSF10C=[0.01;0.03;0.06;0.10;0.20;0.27;t_fit.kd_prot{j};0.5;1];
kp_TNFRSF10C=repmat(t_fit.kp_prot{j}, 1, 9);
kdel_TNFRSF10C=repmat(t_fit.kdel_prot{j}, 1, 9);

for i=1:length(kd_TNFRSF10C)
RNA_fit_puls = fit([-1 0:9]', [1 t_fit.FC_RNA_p{j}(1:10)]','linearinterp');
%[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_TNFRSF10C(i),kd_TNFRSF10C(i),kdel_TNFRSF10C(i)], RNA_fit_puls,  t_fit.FC_Protein_s{j}(1:10), time);
[merr, t_fit_Fig4_protein.c_varykd{i}] = model_protein([kp_TNFRSF10C(i),kd_TNFRSF10C(i),kdel_TNFRSF10C(i)], RNA_fit_puls, kp_TNFRSF10C(i)*ones(length(time),1)/kd_TNFRSF10C(i), time);
t_fit_Fig4_protein.FC_varykd{i}=t_fit_Fig4_protein.c_varykd{i}./repmat(t_fit_Fig4_protein.c_varykd{i}(1),1,length(t_fit_Fig4_protein.c_varykd{i}))
t_fit_Fig4_protein.log2FC_varykd{i}=log2(t_fit_Fig4_protein.FC_varykd{i});
end
fig=figure
set(fig, 'Position',[100 100 1800 300])
for i=1:length(kd_TNFRSF10C)
subplot(1,9,i)
if kd_TNFRSF10C(i)==t_fit.kd_prot{j}
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'linewidth',1,'color','k')
xticks([0:9])
xlim([0,9])
ylim([0,3])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
txt=sprintf('kd fitted=%.2f',kd_TNFRSF10C(i));
text(2,1.7,txt,'FontSize', 9);
else
plot(0:9,t_fit_Fig4_protein.log2FC_varykd{i},'color','b')
xticks([0:9])
xlim([0,9])
ylim([0,3])
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
xlabel('Time (h)')
if i==1
ylabel('TNFRSF10C protein log2(FC)')
end
txt=sprintf('kd=%.2f',kd_TNFRSF10C(i));
text(2,1.7,txt,'FontSize', 9);
end
end
set(gcf,'color','w');
orient(fig,'landscape')
print('-bestfit','plot_vary_kd_TNFRSF10C','-dpdf')
