function [k_opt, x_out, r2] = model_fit(p53_WB, mRNA, time)

k0 = [.1 .2];

k_opt = fmincon(@(x) model(x, p53_WB, mRNA, time), k0, -eye(2), [0;0], ...
    [],[],[],[],[],optimoptions(@fmincon, 'Display','off'));

[merr,x_out] = model(k_opt, p53_WB, mRNA, time);
%  r2 = 1-r2;
r2=corr(ToColumn(mRNA),ToColumn(x_out))^2;

end