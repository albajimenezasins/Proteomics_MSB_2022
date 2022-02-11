
clear all
%load ChipSeq Data
load ./Raw_Data/t_p53ChIP_all_HOMERmotif1.mat
%load TPM counts
Tbl = readtable('./Raw_Data/tximport-tpm.csv');
geneID = Tbl.gene;

%load excell file gene map to translate ensemble id to gene name
[geneMapCCNum,geneMapCCTxt,geneMapCC] = xlsread('./Raw_Data/GRCh38_93_genemap.xlsx','GRCh38.93.genemap');
%extract geneName from geneID
for i=1:length(geneID)
    igeneID=find(strcmp(geneMapCCTxt(:,1),geneID{i}));
    geneNamefound(i)=~isempty(igeneID);
    if (~isempty(igeneID))
        %geneName{i}=geneMapCCTxt(igeneID,3);
        geneName(i)=categorical(geneMapCCTxt(igeneID,3));
    end    
end  

TPM=Tbl(geneNamefound,:)
TPM = TPM(:,3:50);
RNAseqDan.gene = geneName;

%% contaminants and unique genes
%cleam from contaminants, adding 'KRT80' cause not in original excel
contaminants=readtable('./Raw_Data/contaminants.xlsx');
[~,idu] = unique(contaminants(:,3))
contaminants = contaminants(idu,:);
CONT=[contaminants.POLG;'KRT80'];
[Lia,Loc] = ismember(RNAseqDan.gene,CONT);
TPM = TPM(~Lia,:);
RNAseqDan.gene = RNAseqDan.gene(:,~Lia);
%unique genes (maybe isoforms, got rid of them repetitions)
[~,idu] = unique(RNAseqDan.gene);
RNAseqDan.gene= RNAseqDan.gene(:,idu);
TPM = TPM(idu,:);
%% re-order TPM data according to samples
%correspondences also to be found in /n/scratch2/Alba/rnaseq_2020_03_06/200117_A00794_0164_BHHLJ5DRXX/toMerge2.csv
% from mac terminal ssh aj158@o2.hms.harvard.edu then go to folder and cat
% toMerge2.csv to visualize file (use FileZila if need to download
% https://wiki.rc.hms.harvard.edu/display/O2/File+Transfer)
SNAmes=TPM.Properties.VariableNames
samples = readtable('./Raw_Data/samples_summary');
SNAmes=TPM.Properties.VariableNames;
this=samples.Var1';
[Lia,Loc]=ismember(this,SNAmes);
a=TPM.Variables;
RNAseqDan.TPM = a(:,Loc);

RNAseqDan.time=table(samples.Var3,categorical(samples.Var5),samples.Var4,...
    'variablenames', {'time' 'cond' 'replicate'})

%% use ChipSeq Data to add info on which is a p53
p53_targets=unique(categorical(t_p53ChIP_all.GeneSymbol));
[Lia_s,Locb_s] = ismember(RNAseqDan.gene,p53_targets);
RNAseqDan.p53target=Lia_s;

%% pulsatile
% Filter out all genes that have no TPM value at any given time point;
% calculate Fold change and FDR (based on replicates)
p_times = [0:9 12 24 48];
%gene_idx = find(all(RNAseqDan.TPM>=1, 2));
gene_idx=1:length(RNAseqDan.gene)

RNAseqDan_p.TPM = NaN(length(gene_idx), length(p_times));
RNAseqDan_p.pval = NaN(length(gene_idx), length(p_times));
RNAseqDan_p.gene = RNAseqDan.gene(gene_idx);
RNAseqDan_p.p53target = RNAseqDan.p53target(gene_idx);


for iT = 1:length(p_times)
     clear temp1
     clear temp2
     temp1=RNAseqDan.TPM(gene_idx,RNAseqDan.time.time==p_times(iT) & RNAseqDan.time.cond=='p' & RNAseqDan.time.replicate==1);
     temp2=RNAseqDan.TPM(gene_idx,RNAseqDan.time.time==p_times(iT) & RNAseqDan.time.cond=='p' & RNAseqDan.time.replicate==2);
     % modification to keep replicates separate
     RNAseqDan_p.TPM1(:,iT)=temp1(:,1);
     RNAseqDan_p.TPM2(:,iT)=temp2(:,1);
     % mean at 0 should include all timepoints 0, here it is done including sample at 0 from replicate 3
     RNAseqDan_p.TPM(:,iT) = mean(RNAseqDan.TPM(gene_idx, ...
        RNAseqDan.time.time==p_times(iT) & ...
        (RNAseqDan.time.cond=='p' | RNAseqDan.time.time==0) ),2);
    % ttest 0 values (3) against the 2 values selected
    for i = 1:length(gene_idx)
        [~, RNAseqDan_p.pval(i,iT)] = ttest2( ...
            RNAseqDan.TPM(gene_idx(i), RNAseqDan.time.time==0), ...
            RNAseqDan.TPM(gene_idx(i), RNAseqDan.time.time==p_times(iT) & ...
            RNAseqDan.time.cond=='p'));
    end
end

RNAseqDan_p.fdr = NaN(length(gene_idx), length(p_times));
for i = 1:length(gene_idx)
    RNAseqDan_p.fdr(i,2:end) = mafdr(RNAseqDan_p.pval(i,2:end), 'BHFDR', 1);
end

RNAseqDan_p.FC = RNAseqDan_p.TPM./repmat(RNAseqDan_p.TPM(:,1), 1, length(p_times));
RNAseqDan_p.FC1 = RNAseqDan_p.TPM1./repmat(RNAseqDan_p.TPM1(:,1), 1, length(p_times));
RNAseqDan_p.FC2 = RNAseqDan_p.TPM2./repmat(RNAseqDan_p.TPM2(:,1), 1, length(p_times));

RNAseqDan_p.time = p_times;

% added
RNAseqDan_p.log2FC = log2(RNAseqDan_p.FC)
for i = 1:length(RNAseqDan_p.gene)
    RNAseqDan_p.FCs(i,:)=smooth(RNAseqDan_p.log2FC(i,:),4,'moving')';
end

for i = 1:length(RNAseqDan_p.gene)
    repl_matrix=[RNAseqDan_p.FC1(i,1:10); RNAseqDan_p.FC2(i,1:10)]';
    rho=corr(repl_matrix);
    RNAseqDan_p.corr(i)=rho(1,2);
end
%% sustained
%same for sustained condition
s_times = [0:9 12 24 48];

RNAseqDan_s.TPM = NaN(length(gene_idx), length(s_times));
RNAseqDan_s.pval = NaN(length(gene_idx), length(s_times));
RNAseqDan_s.gene = RNAseqDan.gene(gene_idx);
RNAseqDan_s.p53target = RNAseqDan.p53target(gene_idx);

for iT = 1:length(s_times)
     clear temp1
     clear temp2    
     clear temp3
     % times points 0, 1 and 2 from pulsatile copied into sustained 
     if(iT==1|iT==2| iT==3)
     temp1=RNAseqDan.TPM(gene_idx,RNAseqDan.time.time==s_times(iT) & RNAseqDan.time.cond=='p' & RNAseqDan.time.replicate==1);
     temp2=RNAseqDan.TPM(gene_idx,RNAseqDan.time.time==s_times(iT) & RNAseqDan.time.cond=='p' & RNAseqDan.time.replicate==2);
     RNAseqDan_s.TPM1(:,iT)=temp1(:,1);
     RNAseqDan_s.TPM2(:,iT)=temp2(:,1);
     else
     temp1=RNAseqDan.TPM(gene_idx,RNAseqDan.time.time==s_times(iT) & RNAseqDan.time.cond=='s' & RNAseqDan.time.replicate==1);
     temp2=RNAseqDan.TPM(gene_idx,RNAseqDan.time.time==s_times(iT) & RNAseqDan.time.cond=='s' & RNAseqDan.time.replicate==2); 
     temp3=RNAseqDan.TPM(gene_idx,RNAseqDan.time.time==s_times(iT) & RNAseqDan.time.cond=='s' & RNAseqDan.time.replicate==3); 
     RNAseqDan_s.TPM1(:,iT)=temp1(:,1);
     RNAseqDan_s.TPM2(:,iT)=temp2(:,1);   
     % extra sample at time point 6
     if(~isempty(temp3))
       RNAseqDan_s.TPM3=temp3;  
     end
     end
     % for now, sample C6 of replicate 3 is included in the mean
     RNAseqDan_s.TPM(:,iT) = mean(RNAseqDan.TPM(gene_idx, ...
        RNAseqDan.time.time==s_times(iT) & ...
        (RNAseqDan.time.cond=='s' | RNAseqDan.time.time<3)),2);  
    
    for i = 1:length(gene_idx)
        [~, RNAseqDan_s.pval(i,iT)] = ttest2( ...
            RNAseqDan.TPM(gene_idx(i), RNAseqDan.time.time==0), ...
            RNAseqDan.TPM(gene_idx(i), RNAseqDan.time.time==s_times(iT) & ...
            (RNAseqDan.time.cond=='s')));        
    end
end

RNAseqDan_s.fdr = NaN(length(gene_idx), length(s_times));
for i = 1:length(gene_idx)
    RNAseqDan_s.fdr(i,2:end) = mafdr(RNAseqDan_s.pval(i,2:end), 'BHFDR', 1);
end
RNAseqDan_s.FC = RNAseqDan_s.TPM./repmat(RNAseqDan_s.TPM(:,1), 1, length(s_times));


RNAseqDan_s.FC1 = RNAseqDan_s.TPM1./repmat(RNAseqDan_s.TPM1(:,1), 1, length(s_times));
RNAseqDan_s.FC2 = RNAseqDan_s.TPM2./repmat(RNAseqDan_s.TPM2(:,1), 1, length(s_times));

RNAseqDan_s.time = s_times;


% added
RNAseqDan_s.log2FC = log2(RNAseqDan_s.FC)
for i = 1:length(RNAseqDan_s.gene)
    RNAseqDan_s.FCs(i,:)=smooth(RNAseqDan_s.log2FC(i,:),4,'moving')';
end

for i = 1:length(RNAseqDan_s.gene)
    repl_matrix=[RNAseqDan_s.FC1(i,1:10); RNAseqDan_s.FC2(i,1:10)]';
    rho=corr(repl_matrix);
    RNAseqDan_s.corr(i)=rho(1,2);
end



%%
save meanRNAseqDan_all_genes_notTPMfilter.mat RNAseqDan_p RNAseqDan_s 

%% plot safe check genes
% MDM2_expr_p=RNAseqDan_p.TPM(find(RNAseqDan_p.gene=='MDM2'),:)
% MDM2_expr_s=RNAseqDan_s.TPM(find(RNAseqDan_s.gene=='MDM2'),:)
% 
% fig=figure(1);
% set(fig,'PaperSize',[10 9], 'PaperPosition',[0 0 10 9]) 
% subplot (5,3,1)
% [hAx,hLine1,hLine2] = plotyy([0:12]',MDM2_expr_p,[0:12]',MDM2_expr_s)
% set(hAx(1),'ytick',[]);
% linkaxes(hAx, 'xy')
% ylabel(hAx(1),'MDM2 Pulsatile','FontSize', 8) % left y-axis 
% ylabel(hAx(2),'MDM2 Sustained','FontSize', 8) % right y-axis



