function e_shan = entropy_profile(seqs)
% entropy_profile - This function calculates the Shannon entropy within
%                   each position in a multiple alignment
%
% Inputs:
%    seqs - a multiple alighnment represented as a matrix with char
%           ACGT-symbols
%
% Outputs:
%    e_shan - a vector of Shannon entropy values corresponding to each
%             position in the alignment
%
% Other m-files required: uniquemy
%
% Author: Anna A. Igolkina
% email address: igolkinaanna11@gmail.com
% Last revision: 01-Jan-2018

e_shan = [];
for i = 1:size(seqs,2)
    [c, ia] = uniquemy(seqs(:, i));
    ia = ia/sum(ia);
    ia(ia == 0) = [];
    e_shan = [e_shan;-sum(ia .* log(ia))];
end

end









