function entropy_shan = e_shan(seqs)


entropy_shan = [];
for i = 1:size(seqs,2)
    [c, ia] = uniquemy(seqs(:, i));
    ia = ia/sum(ia);
    ia(ia == 0) = [];
    entropy_shan = [entropy_shan;-sum(ia .* log(ia))];
end

end









