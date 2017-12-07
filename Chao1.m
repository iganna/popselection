function schao = Chao1(seqs)

[c, ia] = uniquemy(seqs);

n1 = sum(ia == 1);
n2 = sum(ia == 2);
schao = size(c, 1) + n1*(n1-1) / (2 * (n2+1));

end