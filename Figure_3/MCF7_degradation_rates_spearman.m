
clear all
%load our predicted data through modelling
load fit_results_all.mat
%load Tong et al. 2020 degradation rates
Tbl = readtable('./Raw_Data/Tong_2020_Table_S1.xlsx');

b=Tbl.GeneSymbol
a=t_fit.gene
[Lia,Loc]=ismember(a,b)
Loc=nonzeros(Loc);
t_fit.degradation_Tong_h_1(Lia)=Tbl.DegradationRate_h_1_(Loc);
empties = cellfun('isempty',t_fit.kd_prot);
t_fit.kd_prot(empties) = {0};
kd_prot=cell2mat(t_fit.kd_prot)

a=t_fit.degradation_Tong_h_1(t_fit.degradation_Tong_h_1~=0 & kd_prot~=0)
b=kd_prot(t_fit.degradation_Tong_h_1~=0 & kd_prot~=0)
[rho, pval] = corr(a, b, 'type', 'Spearman')

save fit_results_all_database.mat t_fit

filename = 'modeling_results.xlsx';
writetable(t_fit,filename)