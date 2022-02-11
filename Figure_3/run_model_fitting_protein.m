clear all
load meanProteomics_non_scaled.mat
load fit_results.mat
warning off

% add Protein data
[Lia,Loc] = ismember(t_fit.gene,Proteomics_p.gene);
for i=1:height(t_fit)
    if Lia(i)
        t_fit.FC_Protein_p{i}=Proteomics_p.FC(Loc(i),:);
        t_fit.Protein_corr_p{i}=Proteomics_p.corr(Loc(i));
    end
end
[Lia,Loc] = ismember(t_fit.gene,Proteomics_s.gene);
for i=1:height(t_fit)
    if Lia(i)
        t_fit.FC_Protein_s{i}=Proteomics_s.FC(Loc(i),:);
        t_fit.Protein_corr_s{i}=Proteomics_s.corr(Loc(i));
    end
end

%%
% fit Protein_p to mRNA_p data
time = 0:9;
for i=1:height(t_fit)
    if ~isempty(t_fit.FC_Protein_p{i}) & ~isempty(t_fit.FC_RNA_p{i})

RNA_fit = fit([-1 0:9]', [1 t_fit.FC_RNA_p{i}(1:10)]','linearinterp');
[k_opt, t_fit.fitcurve_prot{i}, t_fit.r2_prot{i},output] = ...
        model_fit_protein(RNA_fit, t_fit.FC_Protein_p{i}(1:10), time);
    
    t_fit.kp_prot{i} = k_opt(1);
    t_fit.kd_prot{i} = k_opt(2);
    t_fit.kdel_prot{i} = k_opt(3);
    end

end
%%
% % fit Protein_s to mRNA_s data
% time = 0:9;
% for i=1:height(t_fit)
%     if ~isempty(t_fit.FC_Protein_s{i}) & ~isempty(t_fit.FC_RNA_s{i})
% 
% RNA_fit = fit([-1 0:9]', [1 t_fit.FC_RNA_s{i}(1:10)]','linearinterp');
% [k_opt, t_fit.FitCurves{i}, t_fit.R2s{i},output] = ...
%         model_fit_protein(RNA_fit, t_fit.FC_Protein_s{i}(1:10), time);
%     
%     t_fit.kp_s{i} = k_opt(1);
%     t_fit.kd_s{i} = k_opt(2);
%     t_fit.kdel_s{i} = k_opt(3);
%     end
% 
% end

%% test on crossing data
% data on sustained and pulsatile (RNA+Prot) is available for 28 genes
time = 0:9;
for i=1:height(t_fit)

      if ~isempty(t_fit.r2_prot{i}) & ~isempty(t_fit.FC_Protein_s{i})

RNA_fit_sus = fit([-1 0:9]', [1 t_fit.FC_RNA_s{i}(1:10)]','linearinterp');
[merr, t_fit.fitcurvepred_prot{i}] = model_protein([t_fit.kp_prot{i},t_fit.kd_prot{i},t_fit.kdel_prot{i}], RNA_fit_sus,  t_fit.FC_Protein_s{i}(1:10), time);
t_fit.r2pred_prot{i}=corr(ToColumn(t_fit.FC_Protein_s{i}(1:10)),ToColumn(t_fit.fitcurvepred_prot{i}))^2;
      end
      
%           if ~isempty(t_fit.R2s{i}) & ~isempty(t_fit.FC_Protein_p{i})
% 
% RNA_fit_puls = fit([-1 0:9]', [1 t_fit.FC_RNA_p{i}(1:10)]','linearinterp');
% [merr, t_fit.fitcurve_prot_froms{i}] = model_protein([t_fit.kp_s{i},t_fit.kd_s{i},t_fit.kdel_s{i}], RNA_fit_puls,  t_fit.FC_Protein_p{i}(1:10), time);
% t_fit.r2_prot_froms{i}=corr(ToColumn(t_fit.FC_Protein_p{i}(1:10)),ToColumn(t_fit.fitcurve_prot_froms{i}))^2;
%           end
end

%% add max of FC proteins
for i=1:height(t_fit)
    a=cell2mat(t_fit.FC_Protein_p(i,:))
    b=cell2mat(t_fit.FC_Protein_s(i,:))
    if ~ isempty(a)
    t_fit.max_FC_Protein_p(i)=max(a(1:10));
    end
    if ~ isempty(b)
    t_fit.max_FC_Protein_s(i)=max(b(1:10));
    end
end
%%
%table with only info of all RNA s, p, Protein s, p
% data on sustained and pulsatile (RNA+Prot) is available for 28 genes

% t_fit_all_info=t_fit;
% loc=cellfun('isempty', t_fit{:,'r2pred_prot'} );
% t_fit_all_info(loc,:)=[];

% % 52 info on proteinp that could be fitted to mRNAp
% loc=cellfun('isempty', t_fit{:,'r2_prot'} );
% sum(loc==0)
% % 49 info on proteinp that could be fitted to mRNAp + protein s that could
% % be crossed checked
% loc=cellfun('isempty', t_fit{:,'r2pred_prot'} );
% sum(loc==0)

% will they be used? 
% for i=1:length(selected_gene)
% t_fitp.FC_RNA {i} = RNAseq_p.FC(selected_gene2 (i),:);
% t_fitp.FC_RNA_max{i}=max(abs(t_fitp.FC_RNA{i}(1:13)));
% t_fitp.FC_RNA_max_all_times{i}=max(abs(t_fitp.FC_RNA{i}));
% end


save fit_results_all.mat t_fit 