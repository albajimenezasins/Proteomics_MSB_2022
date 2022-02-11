
function [merr, x, t_c, x_c] = model(k, p53_WB, mRNA, time)
x = mRNA;
%x = ones(length(time), 1); 
for it = 2:length(time)
    x(it) = x(it-1) + diff(time(it-[1 0]))* ...
        (k(1)*p53_WB(mean(time(it-[1 0]))) - k(2)*x(it-1) );
end

%model_fit(lfit, t_fit_Prot.p(i,1:10), time);
% diff (X)=[X(2)-X(1) X(3)-X(2) ... X(m)-X(m-1)]
% for finding p53 expression RNAs
% load meanRNAseq.mat
% find(RNAseq_p.gene=='TP53')
% [t_c, x_c]=ode45(@(t,x) k(1)*p53_WB(t-.5)-k(2)*x,[time(1) time(end)], mRNA(1));
% 
% x = ToColumn(interp1(t_c, x_c, time, 'linear'));
% mRNA = ToColumn(mRNA);




merr = sum((x-mRNA).^2)/sum((mRNA-mean(mRNA)).^2);

end