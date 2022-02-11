clear all
load clustered_RNAseqDan.mat
cl_RNAseqDan_pulsatile=cl_RNAseqDan_p
clear cl_RNAseqDan_p
clear cl_RNAseqDan_s
load clustered_RNAseqDan_sustained.mat
cl_RNAseqDan_sustained=cl_RNAseqDan_s
clear cl_RNAseqDan_p
clear cl_RNAseqDan_s
load meanRNAseqDan.mat
load meanProteomics_non_scaled.mat


% % filter for differences at early time-points
% Proteins.sumDIVdiffproteinshortterm=nan(length(Proteins.geness),1); 
% for i=1:length(Proteins.geness)     
%     DIV=Proteomics_p.FC(find(Proteomics_p.gene==Proteins.geness(i)),:) ./ Proteomics_s.FC(find(Proteomics_s.gene==Proteins.geness(i)),:)
%     DIVdiff= 1-DIV;
%     Proteins.sumDIVdiffprotein0to3(i)=abs(sum(DIVdiff(1:4)));
% 
% end




% CHOSE GENES that have an upregulated mRNA (either p or s)
% and for which we have info on protein p and protein s and with (either p or s) upregulated
% CORRECTION
% genes that appear induced in S but are noisy
correction={'NDRG1';'TCF12';'LCOR';'RNF19A';'HECTD4';'ELF3';'TOB1';...
    'BCL3';'AREG';'TFAP2C';'IGBP1';'ZMIZ1';'GAB1';'PPFIBP1';'PPM1H';...
    'SASH1';'SEC31A';'MINK1';'CELSR2';'PSD3'}
[Lia_p,Locb_p] = ismember(cl_RNAseqDan_sustained.gene,correction);
cl_RNAseqDan_sustained.cluster(Lia_p)=0
%BCL3 and AREG are cases of single peak cluster in both p53 conditions that 
%lead to differences in expression in protein, they are not representative
%of high mRNA degradation schematics. BCL3 should have been assigned downregulated 
%cluster so we will take it out of both clusters
correction0={'SPG11';'TMTC2';'GRB7';'BCL3';'AREG'}
[Lia_p,Locb_p] = ismember(cl_RNAseqDan_pulsatile.gene,correction0);
cl_RNAseqDan_pulsatile.cluster(Lia_p)=0

f1=(cl_RNAseqDan_pulsatile.cluster==1)|(cl_RNAseqDan_pulsatile.cluster==2)|(cl_RNAseqDan_pulsatile.cluster==3)
genesf1=cl_RNAseqDan_pulsatile.gene(f1)
f2=(cl_RNAseqDan_sustained.cluster==1)|(cl_RNAseqDan_sustained.cluster==2)|(cl_RNAseqDan_sustained.cluster==3)
genesf2=cl_RNAseqDan_sustained.gene(f2)
% unique is an "or"== takes out the repetitions
%if mRNA not present here means and not in another cluster means not
%expressed (flat) or downregulated (cl different than 1, 2, 3)
genes=unique([genesf1;genesf2])
% then select genes for which we have protein info on both sustained and pulsatile
% intersect is an "and"== gives common list of genes
[Lia_p,Locb_p] = ismember(Proteomics_p.gene,genes);
genes_found_p=Proteomics_p.gene(Lia_p);
[Lia_s,Locb_s] = ismember(Proteomics_s.gene,genes);
genes_found_s=Proteomics_s.gene(Lia_s);
genes_we_want=intersect(genes_found_p,genes_found_s)
[Lia_p,Locb_p] = ismember(Proteomics_p.gene,genes_we_want);
%create Proteins
Proteins.geness=Proteomics_p.gene(Lia_p)
Proteins.FCp=Proteomics_p.FC(Lia_p,:)
Proteins.time=Proteomics_p.time;
[Lia_s,Locb_s] = ismember(Proteomics_s.gene,Proteins.geness);
Proteins.FCs=Proteomics_s.FC(Lia_s,:)
% create mRNA
[Lia_p,Locb_p] = ismember(RNAseqDan_p.gene,Proteins.geness);
mRNA.geness=RNAseqDan_p.gene(Lia_p)
mRNA.FCp=RNAseqDan_p.FC(Lia_p,:)
mRNA.time=RNAseqDan_p.time;
[Lia_s,Locb_s] = ismember(RNAseqDan_s.gene,mRNA.geness);
mRNA.FCs=RNAseqDan_s.FC(Lia_s,:)
% cluster 0 means flat
[Lia_p,Locb_p] = ismember(cl_RNAseqDan_pulsatile.gene,mRNA.geness);
Locb_p=nonzeros(Locb_p)
mRNA.clusterp=zeros(length(mRNA.geness),1)
mRNA.clusterp(Locb_p)=cl_RNAseqDan_pulsatile.cluster(Lia_p,:)
[Lia_s,Locb_s] = ismember(cl_RNAseqDan_sustained.gene,mRNA.geness);
Locb_s=nonzeros(Locb_s)
mRNA.clusters=zeros(length(mRNA.geness),1)
mRNA.clusters(Locb_s)=cl_RNAseqDan_sustained.cluster(Lia_s,:)

for i = 1:length(mRNA.geness)
repl_matrix=[mRNA.FCs(i,1:4); mRNA.FCp(i,1:4)]';
rho=corr(repl_matrix);
mRNA.earlycorrsustpuls(i)=rho(1,2);
end



% mRNA and protein corrections
% appear not induced under P but are and belong to cluster 1 (oscillatory)
correction1={'MRI1';'PTPRU';'ATXN3';'ST6GALNAC2';'LIMK2';'ZFP90';'ASCC3';'NR3C1';...
    'STEAP3';'FBXO22';'TMC7';'SFN';'CASK';'RBMS1';'NCEH1';'CMBL';'ISCU';'ITGA3';'BTBD10';'SRA1';'STK17A';'SYT1';'PLXNB1';'MICALL1';'EPHB4';...
    'MGRN1';'RFX7';'NME1';'RIN1';'TSPAN6';'KITLG';'PHPT1';'DSP';'DHRS7';'SLC7A5';'MEGF9';'GBE1';'CHMP3';'TMEM168';'SPPL3';'POLL';'ZSCAN21';'SAC3D1';'FAM84D';'ACYP2';'HS6ST2';'CDC42BPG';'SLC48A1';'HS6ST2';'FAT1'}
[Lia_p,Locb_p] = ismember(mRNA.geness,correction1);
mRNA.clusterp(Lia_p)=1
% appear not induced under P but are and belong to cluster 2 (single peak)
correction2={'CENPP';'HSPA4L';'RGS12';'CCNK';'DNAJB2';'TUBB6';'CASZ1'}
[Lia_p,Locb_p] = ismember(mRNA.geness,correction2);
mRNA.clusterp(Lia_p)=2
% appear not induced under P but are and belong to cluster 3 (rise) 
correction3={'SCRIB';'SEPT9';'PNPO';'CUEDC1';'BAX';'RABGGTA';'FA2H';'PPL';...
    'MOV10';'LMNA';'EPHX1';'AMZ2';'FAM129B';'S100A10';'BAIAP2L1';'CTSL';'ITFG1';'FDXR';'LRPAP1';'PERP';...
    'MAD1L1';'PTPN3';'RPS19';'PDLIM1';'CAMK2G';'CRADD';'KRT19';'PWWP2B';'SPATS2L';'UBE2H'}
[Lia_p,Locb_p] = ismember(mRNA.geness,correction3);
mRNA.clusterp(Lia_p)=3
%do not appear induced in S but are and belond to rising cluster 2
correction4={'COL9A3';'FBXO33'}
[Lia_p,Locb_p] = ismember(mRNA.geness,correction4);
mRNA.clusters(Lia_p)=2

'ITFG1'
%do not appear induced in S but are and belond to rise & decrease cluster1
correction5={'AREG';'OPTN';'PLK2'}
[Lia_p,Locb_p] = ismember(mRNA.geness,correction5);
mRNA.clusters(Lia_p)=1
%appears in an upregulated sustained cluster but we will assign a
%downregulated one such as cluster 5 n=203 or zeros since we selected


SP_mRNA=mRNA.geness((mRNA.clusterp==1|mRNA.clusterp==2|mRNA.clusterp==3)&(mRNA.clusters==1|mRNA.clusters==2|mRNA.clusters==3))
S_mRNA=mRNA.geness(~(mRNA.clusterp==1|mRNA.clusterp==2|mRNA.clusterp==3)&(mRNA.clusters==1|mRNA.clusters==2|mRNA.clusters==3))
P_mRNA=mRNA.geness((mRNA.clusterp==1|mRNA.clusterp==2|mRNA.clusterp==3)& ~(mRNA.clusters==1|mRNA.clusters==2|mRNA.clusters==3))

%compare mRNA sustained to pulsatile expression and store sumDIVdiffmRNA

mRNA.sumDIVdiffmRNA=nan(length(mRNA.geness),1); 
for i=1:length(mRNA.geness)     
    DIV=RNAseqDan_p.FC(find(RNAseqDan_p.gene==mRNA.geness(i)),:) ./ RNAseqDan_s.FC(find(RNAseqDan_s.gene==mRNA.geness(i)),:)
    DIVdiff= 1-DIV;
    mRNA.sumDIVdiffmRNA(i)=sum(DIVdiff(4:11));

end
mRNA.sumDIVdiffmRNA03=(mRNA.sumDIVdiffmRNA>0.3)
% order according to sumDIVdiffmRNA

[DIV_sorted,sort_order]=sort(mRNA.sumDIVdiffmRNA,'descend');
mRNA_ordered = struct( 'FCp', mRNA.FCp(sort_order,:), ...
    'FCs', mRNA.FCs(sort_order,:), ...
    'time', mRNA.time, ...
    'genes', mRNA.geness(sort_order)', ...
    'clusters', mRNA.clusters(sort_order), ...
    'clusterp', mRNA.clusterp(sort_order), ...
    'sumDIVdiffmRNA', mRNA.sumDIVdiffmRNA(sort_order),...
    'sumDIVdiffmRNA03', mRNA.sumDIVdiffmRNA03(sort_order));
mRNA.sumDIVdiffmRNA03(find(mRNA.geness=='ITFG1'))=0;
for i = 1:length(mRNA.geness)
    mRNA_ordered.FC_division(i,:)=mRNA_ordered.FCs(i,:)./mRNA_ordered.FCp(i,:);  
end

%Protein STRUCTURE
% add prot induced information to Proteomics data set
FDR_thres = .05;
FC_thres1 = 1.16;
FC_thres2 = 1.6;
Proteomics_p.protind=(Proteomics_p.corr(:)>0.7)& any(abs(log2(Proteomics_p.FC(:, 2:10)))>log2(FC_thres1),2)| ...
    (isnan(Proteomics_p.corr(:)) & any(abs(log2(Proteomics_p.FC(:, 2:10)))>log2(FC_thres2),2))| ...
    ((Proteomics_p.corr(:)>0.7) & any(abs(log2(Proteomics_p.FC(:, 2:10)))>log2(FC_thres2),2));
Proteomics_s.protind=(Proteomics_s.corr(:)>0.7)& any(abs(log2(Proteomics_s.FC(:, 2:10)))>log2(FC_thres1),2)| ...
    (isnan(Proteomics_s.corr(:)) & any(abs(log2(Proteomics_s.FC(:, 2:10)))>log2(FC_thres2),2))| ...
    ((Proteomics_s.corr(:)>0.7) & any(abs(log2(Proteomics_s.FC(:, 2:10)))>log2(FC_thres2),2));

% protein corrections 
% appear induced under P but are not 
correction0={'MAST4'; 'PNPO'; 'SRA1';'DNAJB2';'S100A11';'CRADD'}
[Lia_p,Locb_p] = ismember(Proteomics_p.gene,correction0);
Proteomics_p.protind(Lia_p)=0
% appear induced under S but are not 
correction0={'NME1';'MARVELD2';'FITM2';'CPSF4'}
[Lia_p,Locb_p] = ismember(Proteomics_s.gene,correction0);
Proteomics_s.protind(Lia_p)=0
% do not appear induced under P but are indeed 
correction1={'TNFRSF10A'; 'EPHB3'; 'SLC12A4';'SCRIB';'PPL';'PTPRU';'DSP'}
[Lia_p,Locb_p] = ismember(Proteomics_p.gene,correction1);
Proteomics_p.protind(Lia_p)=1
% do not appear induced under S but are indeed 
% ST6GALNAC2 case: single replicate in S but similar in expression to
% replicate in P, so I considered it induced
correction2={'PPL';'PTPRU';'STK17A';'RND3';'LRPAP1';'ST6GALNAC2';'DSP'}
[Lia_s,Locb_s] = ismember(Proteomics_s.gene,correction2);
Proteomics_s.protind(Lia_s)=1
% add BCL3 here ;'BCL3'
% 'BLOC1S2' and 'PERP' appear induced under S but are not
correction3={'BLOC1S2';'PERP'}
[Lia_s,Locb_s] = ismember(Proteomics_s.gene,correction3);
Proteomics_s.protind(Lia_s)=0


% make SP, P and S groups 
[Lia_p,Locb_p] = ismember(Proteomics_p.gene,Proteins.geness);
Proteins.ind_P=Proteomics_p.protind(Lia_p);
[Lia_s,Locb_s] = ismember(Proteomics_s.gene,Proteins.geness);
Proteins.ind_S=Proteomics_s.protind(Lia_s);
Proteins.ind_S_not_P=Proteins.ind_S&~Proteins.ind_P
Proteins.ind_P_not_S=Proteins.ind_P&~Proteins.ind_S
Proteins.ind_P_ind_S=Proteins.ind_P&Proteins.ind_S
Proteins.not_P_not_S=~Proteins.ind_P&~Proteins.ind_S

sum(Proteins.ind_P_not_S)
sum(Proteins.ind_S_not_P)
sum(Proteins.ind_P_ind_S)
sum(Proteins.not_P_not_S)

criteria=or(Proteins.ind_S_not_P,Proteins.ind_S)
sum(Proteins.ind_S_not_P|Proteins.ind_P_ind_S)

S_proteins=Proteins.geness(Proteins.ind_S_not_P)
P_proteins=Proteins.geness(Proteins.ind_P_not_S)
SP_proteins=Proteins.geness(Proteins.ind_P_ind_S)


% compare protein sustained to pulsatile expression and store sumDIVdiffmRNA
Proteins.sumDIVdiffprotein=nan(length(Proteins.geness),1); 
for i=1:length(Proteins.geness)     
    DIV=Proteomics_p.FC(find(Proteomics_p.gene==Proteins.geness(i)),:) ./ Proteomics_s.FC(find(Proteomics_s.gene==Proteins.geness(i)),:)
    DIVdiff= 1-DIV;
    Proteins.sumDIVdiffprotein(i)=sum(DIVdiff(4:10));
    Proteins.sumDIVdiffprotein4to7(i)=abs(sum(DIVdiff(5:8)));
    Proteins.sumDIVdiffprotein24(i)=sum(DIVdiff(4:11));

end

Proteins.sumDIVdiffproteinshortterm=nan(length(Proteins.geness),1); 
for i=1:length(Proteins.geness)     
    DIV=Proteomics_p.FC(find(Proteomics_p.gene==Proteins.geness(i)),:) ./ Proteomics_s.FC(find(Proteomics_s.gene==Proteins.geness(i)),:)
    DIVdiff= 1-DIV;
    Proteins.sumDIVdiffprotein0to3(i)=abs(sum(DIVdiff(1:4)));

end


for i=1:length(Proteins.geness)     
Proteins.ratio(i)=Proteins.sumDIVdiffprotein0to3(find(Proteins.geness==Proteins.geness(i)))/Proteins.sumDIVdiffprotein4to7(find(Proteins.geness==Proteins.geness(i)))
end
Proteins.ratio(find(Proteins.geness=='TUBB6'))
Proteins.ratio(find(Proteins.geness=='BLOC1S2'))
Proteins.ratio(find(Proteins.geness=='FDXR'))
Proteins.ratio(find(Proteins.geness=='E2F7'))
sum(Proteins.sumDIVdiffprotein0to3>1)
Proteins.geness(Proteins.sumDIVdiffprotein0to3>1)
Proteins.sumDIVdiffprotein0to3(find(Proteins.geness=='TUBB6')) %1.5
Proteins.sumDIVdiffprotein0to3(find(Proteins.geness=='BLOC1S2')) %1.29
Proteins.sumDIVdiffprotein0to3(find(Proteins.geness=='PERP')) %1.03
Proteins.sumDIVdiffprotein0to3(find(Proteins.geness=='CUEDC1')) %0.34
Proteins.sumDIVdiffprotein0to3(find(Proteins.geness=='HSPG2')) %0.16



Proteins.geness(Proteins.ratio>1)
figure; hist(Proteins.ratio)

Proteins.sumDIVdiffprotein03=(Proteins.sumDIVdiffprotein>0.3)
% order according to sumDIVdifprotein
% winner from 3h to 24h)
%[DIV_sorted,sort_order]=sort(Proteins.sumDIVdiffprotein,'descend');
[DIV_sorted,sort_order]=sort(Proteins.sumDIVdiffprotein24,'descend');
%[DIV_sorted,sort_order]=sort(Proteins.FCs(:,10)./Proteins.FCp(:,10),'descend');
%[DIV_sorted,sort_order]=sort(Proteins.FCs(:,11)./Proteins.FCp(:,11),'descend');


for i = 1:length(Proteins.geness)
repl_matrix=[Proteins.FCs(i,1:4); Proteins.FCp(i,1:4)]';
rho=corr(repl_matrix);
Proteins.earlycorrsustpuls(i)=rho(1,2);
end


Proteins_ordered = struct('ind_P_not_S', Proteins.ind_P_not_S(sort_order,:), ...
    'ind_S_not_P', Proteins.ind_S_not_P(sort_order,:), ...
    'ind_P_ind_S', Proteins.ind_P_ind_S(sort_order,:), ...
    'ind_P', Proteins.ind_P(sort_order,:), ...
    'ind_S', Proteins.ind_S(sort_order,:), ...
    'FCp', Proteins.FCp(sort_order,:), ...
    'FCs', Proteins.FCs(sort_order,:), ...
    'time', Proteins.time, ...
    'genes', Proteins.geness(sort_order)', ...
    'sumDIVdiffprotein', Proteins.sumDIVdiffprotein(sort_order),...
    'sumDIVdiffprotein24', Proteins.sumDIVdiffprotein24(sort_order));


for i = 1:length(Proteins.geness)
    Proteins_ordered.FC_division(i,:)=Proteins_ordered.FCs(i,:)./Proteins_ordered.FCp(i,:);  
end

[Lia,Loc]=ismember(Proteins_ordered.genes,S_proteins)
namesGenes1=cellstr(Proteins_ordered.genes(Lia))
[Lia,Loc]=ismember(Proteins_ordered.genes,P_proteins)
namesGenes2=cellstr(Proteins_ordered.genes(Lia))
[Lia,Loc]=ismember(Proteins_ordered.genes,SP_proteins)
namesGenes3=cellstr(Proteins_ordered.genes(Lia))


% classification mRNA
% S&P
% high kd mRNA belongs to cluster 1 or 2 (transient) mRNAp and clusters 1, 2, 3
criteria= (mRNA.clusterp==1 |mRNA.clusterp==2) & (mRNA.clusters==1 |mRNA.clusters==2|mRNA.clusters==3)
mRNA_SP_highkd=mRNA.geness(criteria)
criteria=mRNA.clusterp==3 & (mRNA.clusters==1 |mRNA.clusters==2|mRNA.clusters==3) & mRNA.sumDIVdiffmRNA03
mRNA_SP_lowkd=mRNA.geness(criteria)
criteria=mRNA.clusterp==3 & (mRNA.clusters==1 |mRNA.clusters==2|mRNA.clusters==3) & ~mRNA.sumDIVdiffmRNA03
mRNA_SP_lowAT=mRNA.geness(criteria)
% P noisy (search_dynamics_mRNA_P)
mRNA_P=P_mRNA
% S (how to distinguish persistence detector vs high and low kd) is a matther of how "flat" the mRNA is
mRNA_S_IFF=S_mRNA

% classification protein
% MDM2 not included in Protein_SP_highkd
Protein_SP_lowkd={'BBC3','EPHA2','TRAF4','RRM2B','TIGAR','IKBIP','ISCU','ITGA3','ARVCF','HSPA4L','BCL2L1','PLCD3','KLF3','EPHB4','PCNA'}
Protein_SP_highkd={'TNFRSF10C','CDKN1A','SESN2','AEN','GDF15','E2F7','FOSL1','PTPRU','TNFRSF10A','PHLDA3','STK17A','CASZ1','RAP2B','PI4K2A','STEAP3','FBXO22','CSNK1G1'}
Protein_SP_lowAT={'POLH','IRF2BPL','PLXNB2','FDXR','GRHL3','SCRIB','PPL','DCP1B','FHL2','PALLD','EPHB3',...
    'XPC','DDB2','RIN1','DGKA','SLC12A4','SLC30A1','SYTL1','AMZ2','SDC4','BTBD10','SYT1','MGRN1','TMEM30A','RFX7','LRPAP1','ST6GALNAC2','ZSCAN21','DSP','FAT1','PTK7'}
Protein_no_class={'PERP','BLOC1S2','CUEDC1','HSPG2','TUBB6'}
% it seemed that 'HS6ST2','MGST2' could be induced in S
% but they have bad correlation between replicates
Protein_S_IFF={'SESN1','NINJ1','PRKAB1','CCNK','MEGF9'}
Protein_S_mediumAT={'MKNK2','HOXC13','KCNN4','LMNA','EPHX1'}
Protein_S_SR={'RAD51C','KDM4B','ARID3A','RHOD'}
Protein_S_SR2={'CTSL','ITFG1','RND3'}

classified_proteins=[Protein_SP_lowkd,Protein_SP_highkd,Protein_SP_lowAT,Protein_S_IFF,Protein_S_mediumAT,Protein_S_SR,Protein_S_SR2,Protein_no_class]
save Fig4_classified_proteins.mat classified_proteins Proteins Proteins_ordered mRNA mRNA_ordered

%checkings
%ismember('CSNK1G1',mRNA_SP_highkd)
%ismember(Protein_SP_lowAT,mRNA_SP_lowkd)
% Proteomics_p.protind(find(Proteomics_p.gene=='CUEDC1'))
% Proteomics_s.protind(find(Proteomics_s.gene=='CUEDC1'))
% Proteomics_s.corr(find(Proteomics_s.gene=='BLOC1S2'))
% Proteomics_p.protind(find(Proteomics_p.gene=='BLOC1S2'))
% Proteomics_s.protind(find(Proteomics_s.gene=='PERP'))
% cl_RNAseqDan_pulsatile.cluster(find(cl_RNAseqDan_pulsatile.gene=='AREG'))
% ismember('BCL3',mRNA_SP_highkd)
%mRNA.clusters(find(mRNA.geness=='FAT1'))
% a=intersect(genesf1,genesf2)
% a=setdiff(b,SP_mRNA)
% s_only=categorical([Protein_S_IFF,Protein_S_mediumAT,Protein_S_SR])
% Protein_highAT=cellstr(Proteins.geness(Proteins.not_P_not_S))
% selected_gene={'PERP','BLOC1S2','CUEDC1','HSPG2','TUBB6'}
% [Lia,Loc]=ismember(selected_gene,Proteomics_s.gene);
% Proteomics_s.corr(Loc)
% Proteomics_s.protind(Loc)'
% [Lia,Loc]=ismember(selected_gene,Proteomics_p.gene);
% Proteomics_p.corr(Loc)
% Proteomics_p.protind(Loc)'
% ismember(selected_gene,mRNA_SP_highkd)
% ismember(selected_gene,mRNA_SP_lowkd)
% ismember('ITFG1',mRNA_SP_lowAT)
% selected_gene=setdiff(categorical(classified_proteins),mRNA.geness)



%TABLE FINAL NUMBERS
%a=17 (MDM2 already out)
a=intersect(mRNA_SP_highkd,Protein_SP_highkd)
%b=13
b=intersect(mRNA_SP_highkd,Protein_SP_lowkd)
%c=23
c=intersect(mRNA_SP_highkd,Protein_SP_lowAT)
%d=5
d=intersect(mRNA_SP_highkd,Protein_S_IFF)
%e=2
e=intersect(mRNA_SP_lowkd,Protein_SP_lowkd)
%f=8
f=intersect(mRNA_SP_lowkd,Protein_SP_lowAT)
%g=5
g=intersect(mRNA_SP_lowkd,Protein_S_mediumAT)
%h=3
h=intersect(mRNA_SP_lowAT,Protein_S_SR2)
%i=4
i=intersect(mRNA_S_IFF,Protein_S_SR)


%find example high AT from mRNA_SP_highkd
list1=mRNA_SP_highkd
list2=[a,b,c,d]
hits=setdiff(list1,list2)
list1=mRNA_SP_lowAT
list2=[h]
hits=setdiff(list1,list2)
list1=mRNA_S_IFF
list2=[i]
hits=setdiff(list1,list2)
