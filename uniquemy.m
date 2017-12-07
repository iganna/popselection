function [c, ia, ic] = uniquemy(seqs)

[c, ~, ic] = unique(seqs, 'rows');
clear ia
for i = 1:size(c,1)
    ia(i,1) = sum(ic == i);
end