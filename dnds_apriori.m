function [n_site, s_site, n_subst, s_subst] = dnds_apriori()
% dnds_apriori - This function calculates values required for dn/ds calculation
%
% Inputs: none
%
% Outputs:
%    n_site - the matrix (64x64) containing at (i,j) position a minimum number 
%    of non-synonymous substitutios which is required to obtain j-th codon
%    from i-th
%    s_site - the matrix (64x64) containing at (i,j) position a minimum number 
%    of synonymous substitutios which is required to obtain j-th codon
%    from i-th
%    n_subst - a vector containing at i-th position
%    a number of single non-synonymous substitutions within i-th codon
%    s_subst - a vector containing at i-th position
%    a number of single synonymous substitutions within i-th codon
%
% Other m-files required: uniquemy, dnds_apriori, change_encoding
%
% Author: Anna A. Igolkina
% email address: igolkinaanna11@gmail.com
% Last revision: 01-Jan-2018


%%
s = [];
for i = 0:63
    j = i;
    s1 = mod(j, 4);
    j = floor(j/4);
    s2 = mod(j, 4);
    j = floor(j/4);
    s3 = mod(j, 4);    
    s = [s;[s1, s2, s3]];
end
s_idx = s+1;
clear s
s_nt = change_encoding(s_idx,1);
s_aa = nt2aamy_nod(s_nt,1);
%% Choose the model
% 'ACGT'
% Jukes-Cantor
ntmx = [0 1/3 1/3 1/3
    1/3 0 1/3 1/3
    1/3 1/3 0 1/3
    1/3 1/3 1/3 0];
ntmx = ntmx*3;

% Kimura Model (transtion - transversion)

a = 2.0192;
b = 1.2938;

ntmx = [0 b b a
    b 0 a b
    b a 0 b
    a b b 0] / (a+b+b); % normalised transitions and transversions


%% Calculate the number of syn/nonsyn sites within each codon
% Here we considet only one-nucleotide-mutation situation.

n_site = zeros(64, 64);
s_site = zeros(64, 64);
for i = 1:64
    if s_aa(i) == '*'
        n_site(i,i) = 3;
        s_site(i,i) = 0;
        continue;
    end
    
    n_by_pos = []; % number of nonsyn when the mutation orrures in one position
    s_by_pos = [];
    for pos = 1:3 % vary each position
        nt1 = s_idx(i, pos); % reference nucleotide
        n_in_pos_by_nt = [];
        s_in_pos_by_nt = [];
        for nt2 = 1:4 % vary each possible nucleotide
            if nt1 == nt2
                continue;
            end
            new_codon_idx = s_idx(i,:);
            new_codon_idx(pos) = nt2;
            new_codon_nt = change_encoding(new_codon_idx, 1);
            

            aa1 = s_aa(i);
            aa2 = nt2aamy_nod(new_codon_nt, 1);
            if aa2 == '*'
                continue;
            end          
            
            if aa1 ~= aa2
                n_in_pos_by_nt = [n_in_pos_by_nt, ntmx(nt1, nt2)];
            else
                s_in_pos_by_nt = [s_in_pos_by_nt, ntmx(nt1, nt2)];            
            end
            clear aa1 aa2
        end
        substitutions_in_pos = sum([n_in_pos_by_nt, s_in_pos_by_nt]);
        if isempty(n_in_pos_by_nt)
            n_in_pos_by_nt = 0;
        end
        if isempty(s_in_pos_by_nt)
            s_in_pos_by_nt = 0;
        end        
        n_by_pos = [n_by_pos, sum(n_in_pos_by_nt) / substitutions_in_pos];
        s_by_pos = [s_by_pos, sum(s_in_pos_by_nt) / substitutions_in_pos];
        clear n_in_pos_by_nt s_in_pos_by_nt
    end
    n_site(i,i) = sum(n_by_pos);
    s_site(i,i) = sum(s_by_pos);
    clear n_by_pos s_by_pos
    clear nt1 nt2
end

 
 

%% Number of substitutions from one codon to another






calculated = zeros(64) ~= 0;
n_subst = zeros(64);
s_subst = zeros(64);

for i = 1:64
    if s_aa(i) == '*'
        continue;
    end    
    for j = 1:64
        if s_aa(j) == '*'
            continue;
        end  
        if i == j
            continue;
        end
        pos = 1:3;
        pos = pos(s_nt(i,:) ~= s_nt(j,:));
        ways = perms(pos);
        n_subst_by_way = [];
        s_subst_by_way = [];        
        n_site_by_way = [];
        s_site_by_way = [];        
        way_with_stop_codon = [];
        for w = 1:size(ways,1)
            way = ways(w, :);
            s_nt_start = s_nt(i,:);
            s_aa_start = nt2aamy_nod(s_nt_start, 1);
            s_nt_end = s_nt(j,:);
            
            for step = 1:length(way)
                s_idx_start = sum((change_encoding(s_nt_start, 0)-1) .* [1 4 16]) + 1;
                n_site_by_way(w, step) = n_site(s_idx_start, s_idx_start);
                s_site_by_way(w, step) = s_site(s_idx_start, s_idx_start);
                
                s_nt_next = s_nt_start;
                s_nt_next(way(step)) = s_nt_end(way(step));
                s_aa_next = nt2aamy_nod(s_nt_next, 1);
                if s_aa_next == '*'
                    way_with_stop_codon = [way_with_stop_codon, w];
                    n_subst_by_way(w,step:length(way)) = nan;
                    s_subst_by_way(w,step:length(way)) = nan;
                    break;
                end
                if s_aa_next == s_aa_start
                    n_subst_by_way(w, step) = 0;
                    s_subst_by_way(w, step) = 1;
                else
                    n_subst_by_way(w, step) = 1;
                    s_subst_by_way(w, step) = 0;
                end
                s_nt_start = s_nt_next;
                s_aa_start = s_aa_next;
            end
        end
        n_subst_by_way(way_with_stop_codon,:) = [];
        s_subst_by_way(way_with_stop_codon,:) = [];
        n_site_by_way(way_with_stop_codon,:) = [];
        s_site_by_way(way_with_stop_codon,:) = [];
        
        n_subst(i, j) = mean(n_subst_by_way(:)) * length(pos);
        s_subst(i, j) = mean(s_subst_by_way(:)) * length(pos);
        n_site(i, j) = mean(n_site_by_way(:));
        s_site(i, j) = mean(s_site_by_way(:));
    end
end



% s_site(s_site == 0) = 1;

