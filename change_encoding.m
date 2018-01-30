function b = change_encoding(a, flag)
% change_encoding - This function change the enciding of nucleatide
% alignment
%
% Inputs: 
%    a - a nuclotice alignment represented by a matrix with ACGT- symbols 
%        or a matrix with 12340 symbols    
%        case of ACGT is not important
%    flag - an instruction what changing should be performed. 
%           if a flag value equals to zero then (ACGT-) -> (12340)
%           if a flag value is not zero (one, for example), then
%                 (12340) -> (ACGT-)
%
% Outputs:
%    a - the matrix with the alignment in another encoding
%
% Other m-files required: none
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