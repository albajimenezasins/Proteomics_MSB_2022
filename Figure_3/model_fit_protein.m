function [k_opt, x_out, r2,output] = model_fit_protein(mRNA, protein_data, time)

k0 = [.1 .2 1];
options=optimoptions('fmincon'); 
options.Algorithm='interior-point';
options.MaxIterations=10000;
options.ConstraintTolerance= 1.0000e-012;
options.OptimalityTolerance= 1.0000e-012;
options.Display='off';


%'interior-point' active-set', 'sqp', 'sqp-legacy', 'trust-region-reflective'
%options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
% introduce upper bound for kdel of 3 hours
lb=[0;0;0];
ub=[Inf;Inf;4];

% [k_opt,a,b,output] = fmincon(@(x) model_protein(x, mRNA, protein_data, time), k0, -eye(3), [0;0;0], ...
%     [],[],[],[],[],options);
[k_opt,a,b,output] = fmincon(@(x) model_protein(x, mRNA, protein_data, time), k0, -eye(3), [0;0;0], ...
    [],[],lb,ub,[],options);
[merr,x_out] = model_protein(k_opt, mRNA, protein_data, time);
r2=corr(ToColumn(protein_data),ToColumn(x_out))^2;
end