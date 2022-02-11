clear all
load clustered_RNAseqDan.mat
%% gene info on clusters 1, 2, 3 of pulsatile
t_fit=table;
t_fit.gene=cl_RNAseqDan_p.gene(cl_RNAseqDan_p.cluster==1|cl_RNAseqDan_p.cluster==2|cl_RNAseqDan_p.cluster==3)
t_fit.cluster=cl_RNAseqDan_p.cluster(cl_RNAseqDan_p.cluster==1|cl_RNAseqDan_p.cluster==2|cl_RNAseqDan_p.cluster==3)
%TPM=cl_RNAseqDan_p.TPM(cl_RNAseqDan_p.cluster==1|cl_RNAseqDan_p.cluster==2|cl_RNAseqDan_p.cluster==3,:)
%TPM_sust=cl_RNAseqDan_s.TPM(cl_RNAseqDan_p.cluster==1|cl_RNAseqDan_p.cluster==2|cl_RNAseqDan_p.cluster==3,:)
FC=cl_RNAseqDan_p.FoldChange(cl_RNAseqDan_p.cluster==1|cl_RNAseqDan_p.cluster==2|cl_RNAseqDan_p.cluster==3,:)
FC_sust=cl_RNAseqDan_s.FoldChange(cl_RNAseqDan_p.cluster==1|cl_RNAseqDan_p.cluster==2|cl_RNAseqDan_p.cluster==3,:)

for a=1:height(t_fit)
%t_fit.TPM{a}=TPM(a,:)
t_fit.FC_RNA_p{a}=FC(a,:)
%t_fit.TPM_sust{a}=TPM_sust(a,:)
t_fit.FC_RNA_s{a}=FC_sust(a,:)
end
%%
%t_fit.TPM=cl_RNAseqDan_p.TPM(cl_RNAseqDan_p.cluster==1|cl_RNAseqDan_p.cluster==2|cl_RNAseqDan_p.cluster==3,:)
%t_fit.TPM_sust=cl_RNAseqDan_s.TPM(cl_RNAseqDan_p.cluster==1|cl_RNAseqDan_p.cluster==2|cl_RNAseqDan_p.cluster==3,:)

% incorporate protein values HERE?

% WB p53 values
p53_puls=[1 4.37 12.78 5.20 1.98 1.50 1.11 3.20 4.01 1.28 1.43 1.85 2.01];
p53_ps=smooth(p53_puls,3);
p53p=[1; p53_ps];
p53_sus=[1 4.37 12.78 15.74 25.69 22.41 21.50 36.88 25.02 37.43 29.85 15.41 22.63];
p53_ss=smooth(p53_sus,3);
p53s=[1; p53_ss];
% actually using the non-smoothed ones
lfit = fit([-1 0:10]', [1 p53_puls(1:11)]','linearinterp');
lfit_sus = fit([-1 0:10]', [1 p53_sus(1:11)]','linearinterp');
% model by gene
time = 0:9;
 for i=1:height(t_fit)
  
    if mod(i,10)==1, disp(i), end
 
    [k_opt, t_fit.fitcurve{i}, t_fit.r2(i)] = ...
        model_fit(lfit, t_fit.FC_RNA_p{i}(time+1), time);

    t_fit.kp(i) = k_opt(1);
    t_fit.kd(i) = k_opt(2);
end

%t_fit = leftjoin(t_fit, t_inter(:,{'gene', ...
    %'hl', 'cluster', 'act_time', 'min_TSS', 'read_count_h0', 'read_count_p0', 'motif_score', 'Max_p53' ,'Max_p53_minus_i'}));

% check previous code run_model_fitting in Tonias folder for the non-peaks



%% test on the sustained data
%t_fit_sus = leftjoin(t_fit(:,{'gene' 'kp' 'kd'}), ...
 %   t_inter(:,{'gene', 'peak_reads_sust' 'TPM_sust'}));
%t_fit_sus =t_fit;

% model by gene
for i=1:height(t_fit)
    if mod(i,10)==1, disp(i), end

    [r2, t_fit.fitcurvepred{i}] = ...
        model([t_fit.kp(i) t_fit.kd(i)], lfit_sus, ...
        t_fit.FC_RNA_s{i}(time+1), time);

t_fit.r2pred(i) = corr( t_fit.FC_RNA_s{i}(time+1)',t_fit.fitcurvepred{i}')^2;
end
%%
save fit_results.mat t_fit

