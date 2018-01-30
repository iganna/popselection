function pchao = pi_chao1(seqs)
% pi_chao1 - This function calculates the Chao1 diversity index
%
% Inputs:
%    seqs - a multiple alighnment represented as a matrix with char
%           ACGT-symbols
%
% Outputs:
%    pchao - value of the Chao1 index
%
% Other m-files required: uniquemy
%
% Author: Anna A. Igolkina
% email address: igolkinaanna11@gmail.com
% Last revision: 01-Jan-2018

[c, ia] = uniquemy(seqs);

n1 = sum(ia == 1);
n2 = sum(ia == 2);
pchao = size(c, 1) + n1*(n1-1) / (2 * (n2+1));

end