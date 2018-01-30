function m_info = mutual_info(seqs)
% mutual_info - This function calculates the mutual (Shannon) information 
%                 between each pair of positions in the multiple alignment
%
% Inputs:
%    seqs - a multiple alighnment represented as a matrix with char
%           ACGT-symbols
%
% Outputs:
%    mutual_info - matrix of mutual information between 
%                  positions in the alignment 
%
% Other m-files required: uniquemy, entropy_profile
%
% Author: Anna A. Igolkina
% email address: igolkinaanna11@gmail.com
% Last revision: 01-Jan-2018

entropy_shan = e_shan(seqs);

m_info = [];
for i = 1:size(seqs,2)
    for j = (i+1):size(seqs,2)

        [~, ia] = uniquemy(seqs(:, [i,j]));
        ia = ia/sum(ia);
        ia(ia == 0) = [];
        
        m_info(i,j) = entropy_shan(i) + entropy_shan(j) - (-sum(ia .* log(ia)));
        m_info(j,i) = entropy_shan(i) + entropy_shan(j) - (-sum(ia .* log(ia)));
    end
end