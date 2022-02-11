
clear all
load fit_results_all.mat
load meanProteomics_non_scaled_R_P.mat
load protein_clusters_pulsatile.mat 

%% kp and kd per mRNA cluster
mRNA_kp1=t_fit.kp((t_fit.cluster==1));
mRNA_kp2=t_fit.kp((t_fit.cluster==2));
mRNA_kp3=t_fit.kp((t_fit.cluster==3));
mRNA_kd1=t_fit.kd((t_fit.cluster==1));
mRNA_kd2=t_fit.kd((t_fit.cluster==2));
mRNA_kd3=t_fit.kd((t_fit.cluster==3));

Prot_kd1=cell2mat(t_fit.kd_prot((t_fit.cluster==1 & ~cellfun(@isempty,t_fit.r2_prot))));
Prot_kd2=cell2mat(t_fit.kd_prot((t_fit.cluster==2 & ~cellfun(@isempty,t_fit.r2_prot))));
Prot_kd3=cell2mat(t_fit.kd_prot((t_fit.cluster==3 & ~cellfun(@isempty,t_fit.r2_prot))));
Prot_kp1=cell2mat(t_fit.kp_prot((t_fit.cluster==1 & ~cellfun(@isempty,t_fit.r2_prot))));
Prot_kp2=cell2mat(t_fit.kp_prot((t_fit.cluster==2 & ~cellfun(@isempty,t_fit.r2_prot))));
Prot_kp3=cell2mat(t_fit.kp_prot((t_fit.cluster==3 & ~cellfun(@isempty,t_fit.r2_prot))));

%% Protein kp and kd per protein cluster
% oscillatory proteins (from 14 oscillatory, 13 are in cl1, 2, 3 of mRNA)
%a=Proteomics_p.gene(Proteomics_p.protein_cluster3_33_bigger==2);
a=cl1_oscillators;
[Lia,Loc]=ismember(a,t_fit.gene);
Prot_kd_cl1_oscill=cell2mat(t_fit.kd_prot(nonzeros(Loc)));
Prot_kp_cl1_oscill=cell2mat(t_fit.kp_prot(nonzeros(Loc)));
%rising proteins (from 15 rising, 15 should be in cl1, 2, and 3)
c=[cl1_rise,cellstr(cl2_rise)',cellstr(cl3_rise)']
%b=Proteomics_p.gene(Proteomics_p.protein_cluster3_33_bigger==3); 
[Lia,Loc]=ismember(c,t_fit.gene);
Prot_kd_cl2_rising=cell2mat(t_fit.kd_prot(nonzeros(Loc)));
Prot_kp_cl2_rising=cell2mat(t_fit.kp_prot(nonzeros(Loc)));
%b=Proteomics_p.gene(Proteomics_p.protein_cluster3_33_bigger==1); 
b=cellstr(cl1_bump)
[Lia,Loc]=ismember(b,t_fit.gene);
Prot_kd_cl3_bump=cell2mat(t_fit.kd_prot(nonzeros(Loc)));
Prot_kp_cl3_bump=cell2mat(t_fit.kp_prot(nonzeros(Loc)));


%% mRNA kd and mRNA kp per mRNA cluster (Tonias observation)
fig=figure
% subplot(1,2,1)
% C = {mRNA_kd1(:); mRNA_kd2(:); mRNA_kd3(:)};  % <--- Just vertically stack all of your groups here
% grp = cell2mat(arrayfun(@(i){i*ones(numel(C{i}),1)},(1:numel(C))')); 
% boxplot(vertcat(C{:}),grp)
% ylim([0 1])
% title('mRNA degradation per cluster','FontSize' ,8)
% ylabel('mRNA kd')
% xlabel('mRNA cluster')
% set(gcf,'color','w');
% subplot(1,2,2)
% C = {mRNA_kp1(:); mRNA_kp2(:); mRNA_kp3(:)};  % <--- Just vertically stack all of your groups here
% grp = cell2mat(arrayfun(@(i){i*ones(numel(C{i}),1)},(1:numel(C))')); 
% boxplot(vertcat(C{:}),grp)
% ylim([0 2])
%ylabel('mRNA kp')
%xlabel('mRNA cluster')
%title('p53 to mRNA production per cluster','FontSize' ,8)
%print(fig,'plot_boxplots_mRNA_kd_kp_per_mRNA_cluster','-dpdf','-r0')
% Protein kd and Protein kp per Protein cluster
%fig=figure
subplot(1,1,1)
C = {Prot_kd_cl1_oscill(:); Prot_kd_cl3_bump(:);Prot_kd_cl2_rising(:)};  % <--- Just vertically stack all of your groups here
grp = cell2mat(arrayfun(@(i){i*ones(numel(C{i}),1)},(1:numel(C))')); 
boxplot(vertcat(C{:}),grp)
ylim([0 0.6])
title('Protein degradation per protein cluster (oscillatory, bump and rising)','FontSize' ,8)
ylabel('Protein kd')
xlabel('Protein cluster')
set(gcf,'color','w');
% subplot(1,2,2)
% C = {Prot_kp_cl1_oscill(:); Prot_kp_cl3_bump(:);Prot_kp_cl2_rising(:)};  % <--- Just vertically stack all of your groups here
% %C = {Prot_kp_cl1_oscill(:); Prot_kp_cl2_rising(:)};% grp = cell2mat(arrayfun(@(i){i*ones(numel(C{i}),1)},(1:numel(C))')); 
% boxplot(vertcat(C{:}),grp)
% ylim([0 0.6])
%ylabel('Protein kp')
%xlabel('Protein cluster')
%title('Protein production per protein cluster (oscillatory vs rising)','FontSize' ,8)
print(fig,'plot_boxplots_mRNA_Protein_kd_kp_FIG_alternative','-dpdf','-r0')