
clear all
load meanRNAseqDan.mat

%Thresholds
FDR_thres = .2;
% correlation filter corr
corr=0.5;
FC_thres = 1.5;
FC_times = [1:9 12];

%Filter genes based on thresholds above
selected_gene = find(any(RNAseqDan_p.fdr<FDR_thres,2) & ...
    any(abs(log2(RNAseqDan_p.FC(:, ismember(RNAseqDan_p.time,FC_times))))>log2(FC_thres),2) &...
   (RNAseqDan_p.p53target'==1)&...
   (RNAseqDan_p.corr>corr)'); 

% correction0={'TRAF4';'PRKAB1'}
% [Lia,Loc] = ismember(RNAseqDan_p.gene,correction0);
% selected_gene=[selected_gene; nonzeros(Loc) ]

%%

sel_RNAseqDan_p = struct('FC', log2(RNAseqDan_p.FC(selected_gene,:)), ...
    'FoldChange', RNAseqDan_p.FC(selected_gene,:), ...    
    'FCs', NaN(length(selected_gene), length(RNAseqDan_p.time)), ...
    'TPM', RNAseqDan_p.TPM(selected_gene,:), ...
    'TPM1', RNAseqDan_p.TPM1(selected_gene,:), ...
    'TPM2', RNAseqDan_p.TPM2(selected_gene,:), ...
    'FC1', RNAseqDan_p.FC1(selected_gene,:), ...
    'FC2', RNAseqDan_p.FC2(selected_gene,:), ...
    'gene', RNAseqDan_p.gene(selected_gene), ...
    'corr', RNAseqDan_p.corr(selected_gene), ...
    'fdr', RNAseqDan_p.fdr(selected_gene,:), ...
    'time', RNAseqDan_p.time);
for i = 1:length(selected_gene)
    sel_RNAseqDan_p.FCs(i,:)=smooth(sel_RNAseqDan_p.FC(i,:),4,'moving')';
end

%including time point 12 for zscore gives a closer result to tonias
Cluster_times = [0:9 12];
z_puls = zscore(sel_RNAseqDan_p.FCs(:,ismember(sel_RNAseqDan_p.time,Cluster_times)),...
    [],2);
%define number of clusters and clustering parameter
n_cluster = 5;
ClusterOpt = [1.3 1000 1e-10 0];

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
[ClusterCenters, ClusterPartition] = fcm(z_puls, n_cluster, ClusterOpt);
[ClusterMembership,ClusterIdx] = max(ClusterPartition);

%% ATTENTION sometimes clusters change, now happens to be cluster 4...
% sum(ClusterIdx==1)
% %32 is single peak (cluster1)
% %38 is increase (cluster2)
% %33 is oscillatory (cluster4)
% selected_gene_cl1=categorical({'CSNK1G1','ZNF33B ','RRM2B','SEMA3E','PIDD1','ARHGEF3','ZNF337',...
%     'TYMSOS','PI4K2A','ANKRA2','DCP1B','ZMAT3','TIGAR','IKBIP','PHLDA3','RIN1','ID3','ACER2','BLOC1S2','CDKN1A','SESN1','TNFRSF10C','TP53INP1'})
% selected_gene_cl3=categorical({'ARVCF','ASTN2'})
% ClusterIdx(ismember(sel_RNAseqDan_p.gene,selected_gene_cl1))=3
% ClusterIdx(ismember(sel_RNAseqDan_p.gene,selected_gene_cl3))=5

%% ATTENTION we are adding a series of important genes whose FC is just below the FC
% max threshold
% should we extend FC to 1.5 and see results?
% best way it to filter by ChipSeq targets, look at FC max distribution,
% then pick a threshold
% adding this genes before clustering does not work so I am adding them
% here
% mRNA corrections these genes appear not induced under P but are and
% should belong to 
% cluster 1 (oscillatory)
% correction1={'MRI1';'PHPT1';'PTPRU';'ATXN3';'ST6GALNAC2';'LIMK2';'ZFP90';'RPS27L';'ASCC3';'NR3C1';...
%     'STEAP3';'FBXO22';'TMC7';'SFN';'DSP';'TRIM32';'IRF2BPL';'CASK';'RBMS1';'ACYP2';'SDC4';'NCEH1';'CMBL';'ISCU';'CASZ1';'ITGA3';'BTBD10';'SRA1';'STK17A'}
% % cluster 2 (single peak)
% correction2={'SYT1';'EPHB3';'CENPP';'DNMBP';'GBE1';'HSPA4L';'CDC42BPG';'TUBB6';'FAT1';...
%     ;'RGS12';'RBM38';'CERS5';'RGS12';'CDC42SE1';'CCNK';'DNAJB2'}
% % cluster 3 (increase)
% correction3={'BAIAP2L1';'SCRIB';'UBE2H';'SEPT9';'PLXNB2';'PNPO';'RASAL1';'PERP';'PTK7';...
%     'LRP5';'DENND2D';'S100A11';'PLXNB1';'KRT19';'BCL2L1';'EPS8L2';'CUEDC1';'BAX';...
%     'MKNK2';'MICALL1';'RABGGTA';'FA2H';'RND3';'TSPAN1';'PPL';'HOXC13';'OCIAD2';'MOV10';'LMNA';'PLCD3';'NTPCR';...
%     ;'CCDC90B';'EPHX1';'AMZ2';'FAM129B';'ARFGAP3';'S100A10'}
% 
% additions=categorical([correction1;correction2;correction3])
% [Lia,Loc] = ismember(additions,RNAseqDan_p.gene);
% selected_gene=[selected_gene; nonzeros(Loc) ]


% 
% sel_RNAseqDan_p = struct('FC', log2(RNAseqDan_p.FC(selected_gene,:)), ...
%     'FoldChange', RNAseqDan_p.FC(selected_gene,:), ...    
%     'FCs', NaN(length(selected_gene), length(RNAseqDan_p.time)), ...
%     'TPM', RNAseqDan_p.TPM(selected_gene,:), ...
%     'TPM1', RNAseqDan_p.TPM1(selected_gene,:), ...
%     'TPM2', RNAseqDan_p.TPM2(selected_gene,:), ...
%     'FC1', RNAseqDan_p.FC1(selected_gene,:), ...
%     'FC2', RNAseqDan_p.FC2(selected_gene,:), ...
%     'gene', RNAseqDan_p.gene(selected_gene), ...
%     'corr', RNAseqDan_p.corr(selected_gene), ...
%     'fdr', RNAseqDan_p.fdr(selected_gene,:), ...
%     'time', RNAseqDan_p.time);
% for i = 1:length(selected_gene)
%     sel_RNAseqDan_p.FCs(i,:)=smooth(sel_RNAseqDan_p.FC(i,:),4,'moving')';
% end
% 
% %add to ClusterIdX
% c1=repmat(4,length(correction1),1)
% c2=repmat(1,length(correction2),1)
% c3=repmat(2,length(correction3),1)
% ClusterIdx=[ClusterIdx,c1',c2',c3']





%% Sort clusters based on activation time (to make the figures look nicer)
clustTime = NaN(1,n_cluster);
act_val= NaN(1,n_cluster);
for i=1:n_cluster 
    %indexes of genes that belong to cluster i
    idx = find(ClusterIdx==i);
    temp=zeros(size(idx,2),1); temp_val=zeros(size(idx,2),1);
    
    % we go through all the genes that belong to cluster i, and for every
    % gene when is the first time that their absolute value goes equal to
    % max temp and their mean value temp_val
    for j=1:size(idx,2)
%     temp(j) = find(abs(puls_fc_s(idx(j),:))>=log2(th),1,'first');
    temp(j) = find(abs(sel_RNAseqDan_p.FC(idx(j),:))>=max(abs(sel_RNAseqDan_p.FC(idx(j),:))),1,'first');
    %temp_val(j)=sel_RNAseqDan_p.FC(idx(j),temp(j));
    temp_val(j)=mean(sel_RNAseqDan_p.FC(idx(j),:));
    
   
    end
    % clustTime average value of activation per cluster, act val intensitie
    % value
    clustTime(i) = mean(temp); act_val(i)=median(temp_val);
    
end

% reorder clusters 1-Nc based on activation time from early to late
[~,orderCluster] = sort((act_val<0)*10+clustTime,'ascend');
ClusterCenters = ClusterCenters(orderCluster,:);
ClusterPartition = ClusterPartition(orderCluster,:);
temp = ClusterIdx;
for i=1:n_cluster 
    ClusterIdx(temp==i) = find(orderCluster==i);
end

ordervalue = 1*ClusterIdx' + sum(abs(sel_RNAseqDan_p.FCs(:,1:13)),2)/(1+max(sum(abs(sel_RNAseqDan_p.FCs(:,1:13)),2)));
% 
[~,sort_order]=sort(ordervalue);
%puls_m_fc_c = sel_RNAseqDan_p.FC(sort_order,:);
%puls_m_fcs_c =  sel_RNAseqDan_p.FCs(sort_order,:);

cl_RNAseqDan_p = struct('FC', sel_RNAseqDan_p.FC(sort_order,:), ...
'FoldChange', sel_RNAseqDan_p.FoldChange(sort_order,:), ...    
'FCs', sel_RNAseqDan_p.FCs(sort_order,:), ...
    'TPM', sel_RNAseqDan_p.TPM(sort_order,:), ...
    'TPM1', sel_RNAseqDan_p.TPM1(sort_order,:), ...
    'TPM2', sel_RNAseqDan_p.TPM2(sort_order,:), ...
    'FC1', sel_RNAseqDan_p.FC1(sort_order,:), ...
    'FC2', sel_RNAseqDan_p.FC2(sort_order,:), ...
    'cluster', ClusterIdx(sort_order)', ...
    'gene', sel_RNAseqDan_p.gene(sort_order)', ...
    'corr', sel_RNAseqDan_p.corr(sort_order)', ...
    'fdr', sel_RNAseqDan_p.fdr(sort_order,:)', ...
    'time', sel_RNAseqDan_p.time, ...
    'ClusterOpt', ClusterOpt);

cl_RNAseqDan_s = struct(...
    'FC', NaN(length(cl_RNAseqDan_p.gene), length(RNAseqDan_s.time)), ...
    'FCs', NaN(length(cl_RNAseqDan_p.gene), length(RNAseqDan_s.time)), ...
    'TPM', NaN(length(cl_RNAseqDan_p.gene), length(RNAseqDan_s.time)), ...
    'cluster', ClusterIdx(sort_order)', ...
    'gene', cl_RNAseqDan_p.gene, ...
    'time', RNAseqDan_s.time);


for i=1:length(cl_RNAseqDan_p.gene)
    cl_RNAseqDan_s.FoldChange(i,:) = RNAseqDan_s.FC(cl_RNAseqDan_p.gene(i)==RNAseqDan_s.gene',:);
    cl_RNAseqDan_s.FC(i,:) = log2(RNAseqDan_s.FC(cl_RNAseqDan_p.gene(i)==RNAseqDan_s.gene',:));
    cl_RNAseqDan_s.FCs(i,:) = smooth(cl_RNAseqDan_s.FC(i,:), 4, 'moving')';
    cl_RNAseqDan_s.TPM(i,:) = RNAseqDan_s.TPM(cl_RNAseqDan_p.gene(i)==RNAseqDan_s.gene,:);
end

%% previously used this now changed before to get good 
% selected_gene_cl1=categorical({'CSNK1G1','ZNF33B ','RRM2B','SEMA3E','PIDD1','ARHGEF3','ZNF337',...
%     'TYMSOS','PI4K2A','ANKRA2','DCP1B','ZMAT3','TIGAR','IKBIP','PHLDA3','RIN1','ID3','ACER2','BLOC1S2','CDKN1A' })
% % create corrected cluster
% cl_RNAseqDan_p.cluster_corrected=cl_RNAseqDan_p.cluster; 
% cl_RNAseqDan_p.cluster_corrected(ismember(cl_RNAseqDan_p.gene,selected_gene_cl1))=1;
%%
save clustered_RNAseqDan.mat cl_RNAseqDan_p cl_RNAseqDan_s 



