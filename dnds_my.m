function [dn, ds, dn_by_pos, ds_by_pos] = dnds_my(seqs, orf, win, type, n_site, s_site, n_subst, s_subst)
% dnds_my - This function calculates an evolutionary statiscics dn/ds
% which is a ratio between relative numers of non-synonymous and synonymous
% substitutions
%
% Inputs:
%    seqs - a sequence alignment presented as a matrix of char
%    orf  - Open reading frame
%    win  - the length of the sliding window. If win is empty or equals to
%           zero the simpson diversity calculates over the whole sequence
%           length
%    type - calculate dn/ds from the "ancestral" (the most abundant
%           haplotype) or not. 
%    n_site, s_site, n_subst, s_subst - values which should be returned
%                                       from dnds_apriori function
%
% Outputs:
%    dn - the number of non-synonymous substitutions per non-synonymous
%    sites
%    ds - the number of synonymous substitutions per synonymous
%    sites
%    dn_by_pos - the number of non-synonymous substitutions per 
%    non-synonymous sites within each codon
%    ds_by_pos - the number of synonymous substitutions per 
%    synonymous sites within each codon
%
% Other m-files required: uniquemy, dnds_apriori
%
% Author: Anna A. Igolkina
% email address: igolkinaanna11@gmail.com
% Last revision: 01-Jan-2018


seqs(:,1:(orf-1)) = [];
seqs = changeCoding(seqs, 0);
clear seqs_aa_idx
for j = 1:3:(size(seqs, 2)-2)
   seqs_codon_idx(:, (j-1)/3 + 1) = ((seqs(:, j:(j+2))-1)*[1; 4; 16]) + 1;
end
seqs_codon_idx(seqs_codon_idx < 0) = 0;
n_codons = size(seqs_codon_idx, 2);
n_seqs = size(seqs_codon_idx, 1);

if isempty(win)
    win = n_codons;   
end
if win == 0
    win = n_codons;
end

% sprintf('window: %i', win)


dn_subst = [];
ds_subst = [];
dn_site = [];
ds_site = [];


% FROM ANCESTRAL NEED TO BE CHECKED
if ~isempty(strfind(type, 'ancestral')) % particial - from ancestral
    [seqs_unique, ia] = uniquemy(seqs_codon_idx);
    seqs_ancestral = seqs_unique(ia == max(ia),:);
    for pos = 1:size(seqs_codon_idx, 2)
        [codons, ia] = uniquemy(seqs_codon_idx(:, pos));
        dn_subst_pos = [];
        ds_subst_pos = [];
        dn_site_pos = [];
        ds_site_pos = [];
        for c = 1:length(codons)
            dn_subst_pos = [dn_subst_pos, n_subst(seqs_ancestral(pos), codons(c)) * ia(c)];
            ds_subst_pos = [ds_subst_pos, s_subst(seqs_ancestral(pos), codons(c)) * ia(c)];
            dn_site_pos = [dn_site_pos, n_site(seqs_ancestral(pos), codons(c)) * ia(c)];
            ds_site_pos = [ds_site_pos, s_site(seqs_ancestral(pos), codons(c)) * ia(c)];
        end
        dn_subst(1,pos) = sum(dn_subst_pos);
        ds_subst(1,pos) = sum(ds_subst_pos);
        dn_site(1,pos) = sum(dn_site_pos);
        ds_site(1,pos) = sum(ds_site_pos);
    end
else
    for pos = 1:size(seqs_codon_idx, 2)
        [codons, ia] = uniquemy(seqs_codon_idx(:, pos));
        dn_subst_pos = [];
        ds_subst_pos = [];
        dn_site_pos = [];
        ds_site_pos = [];


        for c1 = 1:length(codons)
        for c2 = 1:length(codons)

            if c1 == c2 % digonal: no substitutions
                if (codons(c1) + codons(c2)) == 0
                    dn_subst_pos = [dn_subst_pos, 0];
                    ds_subst_pos = [ds_subst_pos, 0];
                    dn_site_pos = [dn_site_pos, 3 * ia(c1) * ia(c1)];
                    ds_site_pos = [ds_site_pos, 0 ];
                else
                    dn_subst_pos = [dn_subst_pos, 0];
                    ds_subst_pos = [ds_subst_pos, 0];
                    dn_site_pos = [dn_site_pos, n_site(codons(c1), codons(c1)) * ...
                        ia(c1) * ia(c1)];
                    ds_site_pos = [ds_site_pos, s_site(codons(c1), codons(c1)) * ...
                        ia(c1) * ia(c1)];
                end
            % supposeing thare is no tolal indel positions
            else
                if (codons(c1) == 0) || (codons(c2) == 0) 
                    codon0 = codons(c1) + codons(c2);
                    dn_subst_pos = [dn_subst_pos,  3 * ia(c1) * ia(c2)];
                    ds_subst_pos = [ds_subst_pos, 0 * ia(c1) * ia(c2)];
                    dn_site_pos = [dn_site_pos, 3 * ia(c1) * ia(c2)];
                    ds_site_pos = [ds_site_pos, 0 * ia(c1) * ia(c2)];
                else
                    dn_subst_pos = [dn_subst_pos, n_subst(codons(c1), codons(c2)) * ...
                          ia(c1) * ia(c2)];
                    ds_subst_pos = [ds_subst_pos, s_subst(codons(c1), codons(c2)) * ...
                        ia(c1) * ia(c2)];
                    dn_site_pos = [dn_site_pos, n_site(codons(c1), codons(c2)) * ...
                        ia(c1) * ia(c2)];
                    ds_site_pos = [ds_site_pos, s_site(codons(c1), codons(c2)) * ...
                        ia(c1) * ia(c2)];
                end
            end
        end
        end
        
        
        dn_subst(1,pos) = sum(dn_subst_pos,2) / n_seqs^2;
        ds_subst(1,pos) = sum(ds_subst_pos,2) / n_seqs^2;
        dn_site(1,pos) = sum(dn_site_pos,2) / n_seqs^2;
        ds_site(1,pos) = sum(ds_site_pos,2) / n_seqs^2;
    end
end


dn = [];
ds = [];
idx = 1:n_codons;
for pos = 1:(n_codons-win+1)
    idx_range = idx((1:win) - 1 + pos);
    dn_up = sum(dn_subst(idx_range));
    ds_up = sum(ds_subst(idx_range));
    dn_down = sum(dn_site(idx_range));
    ds_down = sum(ds_site(idx_range));
    
    if dn_down ~= 0
        dn = [dn, dn_up / dn_down];
    else
        dn = [dn, 0];
    end
    
	if ds_down ~= 0
        ds = [ds, ds_up / ds_down];
    else
        ds = [ds, 0];
	end
    
end


% by position
dn_subst(dn_site == 0) = 0;
dn_site(dn_site == 0) = 1;
ds_subst(ds_site == 0) = 0;
ds_site(ds_site == 0) = 1;

dn_by_pos = dn_subst ./ dn_site;
ds_by_pos = ds_subst ./ ds_site;
end












