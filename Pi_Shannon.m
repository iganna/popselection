function piwin = pi_shannon(seqs, win, mode)
% pi_shannon - This function calculates the Shannon diversity index
%
% Syntax:  piwin = pi_shannon(seqs, win, mode)
%          piwin = pi_shannon(seqs, [], [])
%
% Inputs:
%    seqs - a sequence alignment presented as a matrix of char
%    win  - the length of the sliding window. If win is empty or equals to
%           zero the simpson diversity calculates over the whole sequence
%           length
%    mode - the "shrinkage parameter" paramener to remove sequences with
%           the low abundance.
%           0  : remain all seqs
%           n  : remain seqs with abundance >= n
%           -n : remain seqs with abundance >= n but remove the most abundant
%
% Outputs:
%    piwin - a prifile of Shannon diversity index calculated
%            over a sliding window
%
% Other m-files required: uniquemy
%
% Author: Anna A. Igolkina
% email address: igolkinaanna11@gmail.com
% Last revision: 01-Jan-2018


if isempty(mode)
    mode = 0;
end
if isempty(win)
    win = size(seqs, 2);
end
if win == 0
    win = size(seqs, 2);
end
if isempty(seqs)
    piwin = 0;
    return;
end

clear piwin
for i = 1:(size(seqs, 2) - win + 1)
    [~, ia] = uniquemy(seqs);
    if mode > 0
        ia(ia < mode) = [];
    elseif mode < 0
        ia(ia < -mode) = [];
        ia(ia == max(ia)) = [];
    end
    ia = ia / size(seqs, 1);    
    
    piwin(i) = - sum(ia .* log(ia));
end

end
