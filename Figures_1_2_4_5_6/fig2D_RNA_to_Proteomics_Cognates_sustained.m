clear all
load clustered_RNAseqDan_sustained.mat
load meanProteomics_non_scaled.mat
% 1_check for which mRNAs we have data on Protein
% 2_check for FC protein expression filter 
% 3_cluster proteins independently and plot protein cluster by mRNA clusters 
%% 1  which p53 targets do we have info on proteins n=285 (out of the 410)
% with the 6 clusters in mRNA sustained
[Lia_s,Locb_s] = ismember(cl_RNAseqDan_s.gene,Proteomics_s.gene);
cl_RNAseqDan_s.protinf=Lia_s;
length(cl_RNAseqDan_s.gene(cl_RNAseqDan_s.cluster==1))
length(cl_RNAseqDan_s.gene(cl_RNAseqDan_s.cluster==2))
length(cl_RNAseqDan_s.gene(cl_RNAseqDan_s.cluster==3))
cl1=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protinf&cl_RNAseqDan_s.cluster==1)
cl2=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protinf&cl_RNAseqDan_s.cluster==2)
cl3=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protinf&cl_RNAseqDan_s.cluster==3)



%% 2- which proteins are induced 
% threshold for induced
% a: FC higher 1.1 + corr higher 0.7
% b: FC higher 1.5 + corr Nan or
% c: FC higher 1.5 + corr higher 0.5 
% how many non induced proteins in every cluster?
% we didnt put FDR to be broader
FDR_thres = .05;
FC_thres1 = 1.16;
FC_thres2 = 1.6;

% selected_gene = find((Proteomics_s.corr(:)>0.7)& any(abs(log2(Proteomics_s.FC(:, 2:10)))>log2(FC_thres1),2)| ...
%     (isnan(Proteomics_s.corr(:)) & any(abs(log2(Proteomics_s.FC(:, 2:10)))>log2(FC_thres2),2))| ...
%     ((Proteomics_s.corr(:)>0.5) & any(abs(log2(Proteomics_s.FC(:, 2:10)))>log2(FC_thres2),2)));


protind=(Proteomics_s.corr(:)>0.7)& any(abs(log2(Proteomics_s.FC(:, 2:10)))>log2(FC_thres1),2)| ...
    (isnan(Proteomics_s.corr(:)) & any(abs(log2(Proteomics_s.FC(:, 2:10)))>log2(FC_thres2),2))| ...
    ((Proteomics_s.corr(:)>0.7) & any(abs(log2(Proteomics_s.FC(:, 2:10)))>log2(FC_thres2),2));

%[Lia_s,Locb_s] = ismember(cl_RNAseqDan_s.gene,Proteomics_s.gene);
[Lia_s,Locb_s] = ismember(Proteomics_s.gene,cl_RNAseqDan_s.gene);
id=nonzeros(Locb_s);
cl_RNAseqDan_s.protind=zeros(length(cl_RNAseqDan_s.gene),1)
cl_RNAseqDan_s.protind(id)=protind(Lia_s);

%of induction
%25 induced out of 34 info cl1 % 73%
%45 induced out of 116 info cl2 % (43) 38%
%6 induced out of 62 info cl3% (6) 10%
% total induction 25+45+6 (76) out of 34+116+62 (212) 65%


cl1=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protinf&cl_RNAseqDan_s.protind&cl_RNAseqDan_s.cluster==1)
cl2=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protinf&cl_RNAseqDan_s.protind&cl_RNAseqDan_s.cluster==2)
cl3=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protinf&cl_RNAseqDan_s.protind&cl_RNAseqDan_s.cluster==3)
cl4=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protinf&cl_RNAseqDan_s.protind&cl_RNAseqDan_s.cluster==4)
cl5=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protinf&cl_RNAseqDan_s.protind&cl_RNAseqDan_s.cluster==5)
cl6=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protinf&cl_RNAseqDan_s.protind&cl_RNAseqDan_s.cluster==6)


cl1_o=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protinf&cl_RNAseqDan_s.cluster==1)
cl2_o=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protinf&cl_RNAseqDan_s.cluster==2)
cl3_o=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protinf&cl_RNAseqDan_s.cluster==3)
%TOTAL INDUCTION 34+117+65=216 INFO / 76 INDUCED= 35%
%% 3- cluster prots independently 
% although we are interested in the proteins from cl1, 2, 3, that are
% not in the decaying mRNA clusters
% clustering on all induced proteins 82 doesnt make a difference so we go for that 
% OR 
cl=[cl1;cl2;cl3;cl4;cl5;cl6]
[Lia,Loc] = ismember(Proteomics_s.gene,cl);
clusterIDorder=Proteomics_s.gene(Lia)

z_puls = zscore(Proteomics_s.smoothlog2FC(Lia,1:10),[],2);
z_puls_all_prots = zscore(Proteomics_s.smoothlog2FC(:,1:10),[],2);

n_cluster = 2;
ClusterOpt = [1.3 1000 1e-10 0];

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
[ClusterCenters, ClusterPartition] = fcm(z_puls, n_cluster, ClusterOpt);
[ClusterMembership,ClusterIdx] = max(ClusterPartition);
Proteomics_s.protein_cluster2_82=ClusterIdx;
cl_RNAseqDan_s.protein_cluster2_82=zeros(length(cl_RNAseqDan_s.cluster),1)
[Lia2,Loc2] = ismember(Proteomics_s.gene(Lia),cl_RNAseqDan_s.gene)
cl_RNAseqDan_s.protein_cluster2_82(Loc2)=ClusterIdx;


%% plot protein clusters
n_cluster = max(Proteomics_s.protein_cluster2_82);
colors = [[.69:.15:1 ; 0:.25:.6; 0:.2:.45]'
    [0:.2:.45; 0:.25:.6; .69:.15:1]'];
fig=figure;
set( fig,'PaperSize',[10 3], 'PaperPosition',[0 0 10 3]) 
%
[a,b]=hist(Proteomics_s.protein_cluster2_82,unique(Proteomics_s.protein_cluster2_82))
for i=1:n_cluster
    subplot(1,n_cluster,i)
    plot(0:9,z_puls(Proteomics_s.protein_cluster2_82==i,:),'color', colors(i,:),'linewidth',0.5);
    hold on
    plot(0:9,mean(z_puls(Proteomics_s.protein_cluster2_82==i,:)),'k','linewidth',2)   
    text(2,2,sprintf('cluster %.d: n=%.d',i,a(i)),'FontSize',7)
    xlim([0 9])
    ylim([-3 3])
    set(gca,'xtick',0:9,'xticklabel',0:9,'FontSize',7)   
end
xlabel('time[h]')
[ax,h1]=suplabel('time[h]','x');
[ax,h2]=suplabel('zscore','y'); 
print('fig2_Proteomics_Cognates_clusters_sustained','-dpdf','-r0')
%% plot proteins according to mRNA targets cluster and protein cluster
z_puls = zscore(Proteomics_s.smoothlog2FC(:,1:10),[],2);
fig=figure('Renderer', 'painters', 'Position', [1 400 500 950 ])

for a=1:max(cl_RNAseqDan_s.cluster)
    %for a=1:3
for i=1:max(cl_RNAseqDan_s.protein_cluster2_82)
    subplot(max(cl_RNAseqDan_s.cluster),max(cl_RNAseqDan_s.protein_cluster2_82),i+(a-1)*2)
    genes=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protein_cluster2_82==i&cl_RNAseqDan_s.cluster==a)
    if length(genes)~=0
    [Lia,Loc] = ismember(genes,Proteomics_s.gene);
    Loc=nonzeros(Loc)
    %plot(0:9,Proteomics_s.FC(Loc,[1:10]),'linewidth',0.5);
    plot(0:9,z_puls(Loc,[1:10])','linewidth',0.5);
    ylim([-2.5 2.5])
    hold on
    plot(0:9,mean(z_puls(Loc,[1:10]))','k','linewidth',1)   
    text(5,1,sprintf('n=%.2f',length(Loc)),'FontSize',7)
    title(sprintf('prot cl.%d from mRNA cl.%d',i,a),'FontSize',7) 
    if mod(i+(a-1)*2,2)
        ylabel('zscore')
    end
    xlim([0 9])
    set(gca,'xtick',0:9,'xticklabel',0:9)
    if (i==2 & a<4)     
    x=0:9;
    y=mean(z_puls(Loc,[1:10]))
    P = polyfit(x,y,1);
    slope = P(1)
    text(5,-1,sprintf('slope=%.2f',slope),'FontSize',7)    
    ArrCell = {}
    for c=1:length(Loc)
        y2=z_puls(Loc(c),[1:10]);
        Pol= polyfit(x,y2,1);
        ArrCell= [ { Pol(1) }, ArrCell];
    end
    Protein_slopes{a}=ArrCell    
    end 
    end
   
end
end
xlabel('time [h]')
suplabel('protein expression by mRNA cluster') % this to be on the right side

set(fig,'Units','Inches');
set(gcf,'color','w');
set(fig,'PaperOrientation','portrait');
print(fig,'fig2_RNA_Proteomics_Cognates_sustained','-dpdf','-r0')


%rise proteins
rise_genes=cl_RNAseqDan_s.gene(cl_RNAseqDan_s.protein_cluster2_82==2)
% store cluster info in Proteomics_p.protein_cluster
Proteomics_s.protein_cluster=zeros(length(Proteomics_s.gene),1)
[Lia, Loc]=ismember(rise_genes,Proteomics_s.gene)
Proteomics_s.protein_cluster(Loc)=3;





save clustered_RNAseqDan_R_P_sustained.mat cl_RNAseqDan_p cl_RNAseqDan_s
save meanProteomics_non_scaled_R_P_sustained.mat Proteomics_s
