function a = nt2aamy_nod(nts, ORF)
% nt2aamy_nod - converstor of an alignment of nucleotide sequences to the
%           alignment of amino acid sequences
%
% Inputs:  
%          nts - a multiple alighnment represented as a matrix with char
%                ACGT-symbols
%          ORF - Open Reading Frame: start of the codon (1, 2 or 3)
%
% Outputs:
%    a - amino acid alignment
%
% Other m-files required: none
%
% Author: Anna A. Igolkina
% email address: igolkinaanna11@gmail.com
% Last revision: 01-Jan-2018


nts = upper(nts);
code = ['KAAA'
'NAAC'
'KAAG'
'NAAT'
'TACA'
'TACC'
'TACG'
'TACT'
'RAGA'
'SAGC'
'RAGG'
'SAGT'
'IATA'
'IATC'
'MATG'
'IATT'
'QCAA'
'HCAC'
'QCAG'
'HCAT'
'PCCA'
'PCCC'
'PCCG'
'PCCT'
'RCGA'
'RCGC'
'RCGG'
'RCGT'
'LCTA'
'LCTC'
'LCTG'
'LCTT'
'EGAA'
'DGAC'
'EGAG'
'DGAT'
'AGCA'
'AGCC'
'AGCG'
'AGCT'
'GGGA'
'GGGC'
'GGGG'
'GGGT'
'VGTA'
'VGTC'
'VGTG'
'VGTT'
'*TAA'
'YTAC'
'*TAG'
'YTAT'
'STCA'
'STCC'
'STCG'
'STCT'
'*TGA'
'CTGC'
'WTGG'
'CTGT'
'LTTA'
'FTTC'
'LTTG'
'FTTT'
'----'
];

a = '';
for i = 1:size(nts,1)
     for j = ORF:3:(size(nts, 2))
        if (j+2) > (size(nts, 2))
            continue;
        end
        nt = nts(i,(j+[0:2]));
%         nt = upper(nt);
        ident = [code(:,2) == nt(1), code(:,3) == nt(2), code(:,4) == nt(3)];
        ident = sum(ident, 2) == 3;
        if sum(ident) == 0
            a(i, (j-ORF)/3+1) = '-';
        else
            a(i, (j-ORF)/3+1) = code(ident, 1);    
        end
    end
end

