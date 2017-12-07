function mutual_info = mut_shan(seqs)

entropy_shan = e_shan(seqs);

mutual_info = [];
for i = 1:size(seqs,2)
    for j = (i+1):size(seqs,2)

        [c, ia] = uniquemy(seqs(:, [i,j]));
        ia = ia/sum(ia);
        ia(ia == 0) = [];
        
        mutual_info(i,j) = entropy_shan(i) + entropy_shan(j) - (-sum(ia .* log(ia)));
        mutual_info(j,i) = entropy_shan(i) + entropy_shan(j) - (-sum(ia .* log(ia)));
    end
end