function b = change_encoding(a, flag)
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
if flag == 0
    b = [];
    for i = 1:size(a,1)
        for j = 1:size(a, 2)
            if (a(i,j) == 'A') || (a(i,j) == 'a')
                b(i,j) = 1;
            elseif (a(i,j) == 'C') || (a(i,j) == 'c')
                b(i,j) = 2;
            elseif (a(i,j) == 'G') || (a(i,j) == 'g')
                b(i,j) = 3;                
            elseif (a(i,j) == 'T') || (a(i,j) == 't')
                b(i,j) = 4;                
            elseif (a(i,j) == '-')
                b(i,j) = 0;                
            else
                b(i,j) = -1;                
            end
        end
    end
else
    b = '';
    for i = 1:size(a,1)
        for j = 1:size(a, 2)
            if (a(i,j) == '1') || (a(i,j) == 1)
                b(i,j) = 'A';
            elseif (a(i,j) == '2') || (a(i,j) == 2)
                b(i,j) = 'C';
            elseif (a(i,j) == '3') || (a(i,j) == 3)
                b(i,j) = 'G';                
            elseif (a(i,j) == '4') || (a(i,j) == 4)
                b(i,j) = 'T';                
            elseif (a(i,j) == 0)
                b(i,j) = '-';  
            else
                b(i,j) = 'N';  
            end
        end    
    end


end

end