function [pvalues_modchisquare] = CalculateChisquarePvalues(observed_chisquare,observed_count,df,pvalr)
%CalculateChisquarePvalues calculates the pvalues for a set of chisquare statistics. 
%
%USAGE: pvalues_modchisquare = CalculateChisquarePvalues(observed_chisquare,observed_count,df,pvalr);
%
%INPUTS:
%   observed_chisquare -- a NxN matrix of observed chisquare statistics
%
%   observed_count -- a NxNx4 matrix of observed and expected counts. The
%   second sheet should be observed significant values, and the fourth sheet
%   should be expected significant values. If you are not using
%   significance tests, simply set the second sheet values to be greater
%   than the fourth sheet values.
%
%   df -- the degrees of freedom, usually 1, but can be more if this is 
%   being used outside of CountSignificantEffectsByModules.
%
%   pvalr -- the value used to calculate the false discovery rate (FDR).
%
%OUTPUTS:
%
%   pvalues_modchisquare -- a NxNx2 matrix where the first sheet contains
%   the p values and the second sheet is a matrix of ones,zeros, and
%   negative ones (for significantly greater than chance (1), significantly
%   less than chance (-1), or no significance (0)).
%
%NOTES:
%  
%   This function is called in CountSignificantEffectsByModules
%
%SEE ALSO: CountSignificantEffectsByModules

%VERSION HISTORY%
%
%Version 1.0, initialized by Eric Feczko on 4/30/14

pvalues_modchisquare = zeros(size(observed_chisquare,1),size(observed_chisquare,2),2);
pvalues_modchisquare(:,:,1) = chi2cdf(observed_chisquare,df);

for i = 1:size(pvalues_modchisquare,1)
    for j = 1:size(pvalues_modchisquare,2)
    if pvalues_modchisquare(i,j,1) == 1
        test = 1;
        stringval = '0.';
        while test == 1
            stringval = strcat(stringval,'9');
            chivalue = chi2inv(str2num(stringval),df);
            if chivalue <= observed_chisquare(i,j)
                newpvalue = 1 - str2num(stringval);
            else
                test = 0;
            end
        end
        pvalues_modchisquare(i,j,1) = newpvalue;
    else
        pvalues_modchisquare(i,j,1) = 1 - pvalues_modchisquare(i,j,1);
    end
    end
end
pvalues_list = zeros(size(pvalues_modchisquare,1)^2,1);
pvalues_list(:) = pvalues_modchisquare(:,:,1);
pvalues_list_sorted = sort(pvalues_list);
V = length(pvalues_list_sorted);
I = (1:V)';

cVID = 1;
cVN = sum(1./(1:V));

pID = pvalues_list_sorted(find(pvalues_list_sorted<=I/V*pvalr/cVID, 1, 'last' ))
if(size(pID,1) == 0)
    pID = 0;
    'no connections survive multiple comparison correction'
end
pN = pvalues_list_sorted(find(pvalues_list_sorted<=I/V*pvalr/cVN, 1, 'last' ))
if(size(pN,1) == 0)
    pN =0;
    'no connections survive multiple comparison correction'
end
pvalues_dist = zeros(size(observed_chisquare,1),size(observed_chisquare,2));
for i = 1:size(observed_chisquare,1)
    for j = 1:size(observed_chisquare,2);
        if observed_count(i,j,2) > observed_count(i,j,4)
            pvalues_dist(i,j) = 1;
        else
            pvalues_dist(i,j) = -1;
        end
    end
end
for i = 1:size(observed_chisquare,1)
    for j = 1:size(observed_chisquare,2);
        if pvalues_modchisquare(i,j,1) <= pvalr
            pvalues_modchisquare(i,j,2) = 1*pvalues_dist(i,j);
        end
    end
end

end

