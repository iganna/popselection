function [piwin] = Pi_Shannon(seqs, win, mode)



% 0 : All seqs
% n : more than n
% -n : more than n and withoun the most abundant

if isempty(win)
    win = 1;
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
    
    [c, ia] = uniquemy(seqs);

    
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