function [c, ia, ic] = uniquemy(seqs)
% uniquemy - This function returns the unique sequences and their frequences
%            However, this function can digest also vertical vectors of
%            cells
%
% Syntax:  [c, ia, ic] = uniquemy(seqs)
%
% Inputs:
%    seqs - It can be (1) a sequence alignment presented as a matrix of char
%                     (2) a vertical vector of cells
%
% Outputs:
%    c - unique sequences
%    ia - a vertical vector with frequencies
%    ic - a vector of indexes: seqs = c(ic,:)
%
% Other m-files required: none
%
%
% Author: Anna A. Igolkina
% email address: igolkinaanna11@gmail.com
% Last revision: 01-Jan-2018

%------------- BEGIN CODE --------------

[c, ~, ic] = unique(seqs, 'rows');
clear ia
for i = 1:size(c,1)
    ia(i,1) = sum(ic == i);
end
%------------- END OF CODE --------------

