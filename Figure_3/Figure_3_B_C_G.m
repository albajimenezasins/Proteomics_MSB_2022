clear all
load fit_results_all.mat
load clustered_RNAseqDan_R_P.mat

%% USE PLOT GENERAL AND CALL FIG 2D AND 2E
%% genes with bad R2 despite good FC , Fig2D
%selected_genes={'ZNF219','TNFAIP8', 'BLOC1S2'}
%% genes with bad R2 sustained despite good R2 pulsatile, Fig 2E
%selected_genes={'HSPG2','GSN','MAST4','CHSY1','SIPA1L3'}
%% data and treat zeros first
empties = cellfun('isempty',t_fit.r2_prot);
t_fit.r2_prot(empties) = {0};
empties = cellfun('isempty',t_fit.r2pred_prot);
t_fit.r2pred_prot(empties) = {0};
a=cell2mat(t_fit.r2_prot);
b=t_fit.max_FC_Protein_p;
c=cell2mat(t_fit.r2pred_prot);

%mRNA
a_mRNA=t_fit.r2
c_mRNA=t_fit.r2pred;


%% replace zeros with NaNs
a_corr=a
a_corr(a_corr==0)=NaN
b_corr=b
b_corr(b_corr==0)=NaN
c_corr=c
c_corr(c_corr==0)=NaN
% get rid of NaNs is order to do Spearman correlations
a_c=a_corr
b_c=b_corr
a_c(isnan(a_corr) | isnan(b_corr))=[];
b_c(isnan(a_corr) | isnan(b_corr))=[];
a_d=a_corr
c_d=c_corr
a_d(isnan(a_corr) | isnan(c_corr))=[];
c_d(isnan(a_corr) | isnan(c_corr))=[];
[rho, pval] = corr(a_c, b_c, 'type', 'Spearman')
[rho_r, pval_r] = corr(a_d, c_d, 'type', 'Spearman')

%% sort order
% THIS WAS BETTER DONE BY JUST PLOTTING THE SCATTER PLOT
% % order then on sustained bad
% empties = cellfun('isempty',t_fit.r2_prot);
% t_fit.r2_prot(empties) = {0};
% empties = cellfun('isempty',t_fit.r2pred_prot);
% %t_fit.r2pred_prot(empties) = {NaN};
% t_fit.r2pred_prot(empties) = {10};
% % order first on pulsatile good
% [a,sort_order1]=sort(cell2mat(t_fit.r2_prot),'descend');
% %t_fit_r2_prot_order = t_fit(sort_order,:);
% % order then on sustained bad
% [b,sort_order2]=sort(cell2mat(t_fit.r2pred_prot),'ascend');
% [~,sort]=sort(sort_order1+sort_order2); 
% t_fit_sorted= t_fit(sort,:);
% % genes HSPG2, GSN, MAST4, CHSY1 have good R2 pulsatile and bad R2


%% plot R2 the same 28 proteins than Fig2 then reuse previous plots
fig2_genes=[]
for a=1:3
for i=1:max(cl_RNAseqDan_p.protein_cluster3_33)
    genes=cl_RNAseqDan_p.gene(cl_RNAseqDan_p.protein_cluster3_33==i&cl_RNAseqDan_p.cluster==a)
    fig2_genes=[fig2_genes;genes]
end
end
fig2_genes=fig2_genes(~ismember(fig2_genes,'MAST4'))
%find data
[Lia,Loc]=ismember(t_fit.gene,fig2_genes);
Loc=nonzeros(Loc);
aa=median(t_fit.r2(:))
a=cell2mat(t_fit.r2_prot(Lia));
b=t_fit.max_FC_Protein_p(Lia);
empties = cellfun('isempty',t_fit.r2pred_prot);
t_fit.r2pred_prot(empties) = {0};
c=cell2mat(t_fit.r2pred_prot(Lia));
[rho, pval] = corr(a, b, 'type', 'Spearman')
[rho_r, pval_r] = corr(a, c, 'type', 'Spearman')

%% histogram for Fig
figure
histogram(c,20)
xlabel('R2 Prediction ')
ylabel('number of genes')
set(gcf,'color','w');
ylim([0,16])
set(gca,'fontsize',18)
print('-bestfit','plot_R2_histogram','-dpdf')
%%
figure
subplot(1,2,1)
histogram(a_mRNA,40,'Normalization','probability')
xlabel('R2 mRNA fit ')
ylabel('number of genes')
set(gcf,'color','w');
ylim([0,0.12])
%set(gca,'fontsize',18)
%print('-bestfit','plot_R2_mRNA_histogram','-dpdf')
subplot(1,2,2)
histogram(c_mRNA,30,'Normalization','probability')
xlabel('R2 mRNA pred ')
ylabel('number of genes')
set(gcf,'color','w');
ylim([0,0.22])
%set(gca,'fontsize',18)
print('-bestfit','plot_R2_mRNA_histograms_new_new','-dpdf')
