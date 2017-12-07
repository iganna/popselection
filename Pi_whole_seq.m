function [pi_abs, pi_norm] = Pi_whole_seq(s1, s2)
% pi - absulite value of differences
% normalised = probability
len = size(s1, 2);

pi_abs = 0;
if isempty(s2)

    [s1, cnt1, ~] = uniquemy(s1);
    n1 = length(cnt1); % unique
    
    for j = 1:n1
        for k = (j + 1):n1
            pi_abs = pi_abs + cnt1(j) * cnt1(k) * sum(s1(j,:) ~= s1(k,:));
        end
    end
    
    n1 = sum(cnt1); % total
    pi_norm = pi_abs / (n1 * (n1 - 1) / 2);
    
else

    [s1, cnt1, ~] = uniquemy(s1);
    [s2, cnt2, ~] = uniquemy(s2);
    n1 = length(cnt1); % unique
    n2 = length(cnt2); % unique
    
    for j = 1:n1
        for k = 1:n2
            pi_abs = pi_abs + cnt1(j) * cnt2(k) * sum(s1(j,:) ~= s2(k,:));
        end
    end       
    
    n1 = sum(cnt1); % unique
    n2 = sum(cnt2); % unique
    pi_norm = pi_abs / (n1 * n2);
end

pi_abs = pi_abs / len;
pi_norm = pi_norm / len;


end