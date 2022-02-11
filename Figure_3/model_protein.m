
function [merr, x, t_c, x_c] = model_protein(k, mRNA, protein, time)
%x(1) = protein(1);
x = protein;
for it = 2:length(time)
    x(it) = x(it-1) + k(1)*mRNA(time(it)-k(3)) - k(2)*x(it-1);
    % used to be k(1)*mRNA(mean(time(it-[1 0]))) where mRNA(Y) and Y 0.5 to 8.5 so that timepoint 2 relates to 0.5
end

merr = sum((x-protein).^2)/sum((protein-mean(protein)).^2);

end