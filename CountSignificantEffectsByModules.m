function [module_mat_chisquare module_mat_count module_mat_ratio module_mat_pvalue] = CountSignificantEffectsByModules(m,modules,df,pvalr)
%CountSignificantEffectsByModules will take a binarized matrix of signifiance
%tests (1 - significant, 0 - not significant) and determine whether those
%tests occur within/between specific communities (modules).
%
%USAGE: [module_mat_chisquare, module_mat_count, module_mat_ratio,module_mat_pvalue] = CountSignificantEffectsByModules(m,modules,df,pvalr);
%
%INPUTS:
%
%   m -- A binarized matrix for signifiant (1) and not significant (0)
%   events. Each row/column in the matrix reflects a single region of interest (ROI). 
%   Typically, a binarized matrix results from thresholding a
%   set of significance tests (e.g. t-test or pearson's correlation)
%
%   modules -- a 1-D matrix where the length of the matrix is the same size
%   as the number of ROIs in the binarized matrix. Each item in the matrix
%   is a number that represents the community assignment for the given ROI.
%   The order of items should match the order of the binarized matrix.
%
%   df -- the degrees of freedom for performing the chi-square test. Since
%   the chi-square test here is a 2x1 (significant/non-significant) the
%   degrees of freedom should be 1.
%
%   pvalr -- The threshold for calculating the false-discovery rate for a
%   test. Typically this is set to 0.05/N, where N is the number of
%   matrices that are being clustered.
%
%
%OUTPUTS:
%
%   module_mat_chisquare -- a community X community matrix containing the 
%   chisquare values for the modules.
%
%   module_mat_count -- a community X community X 4 matrix containing the
%   raw counts. The first sheet is the observed number of non-significant effects per
%   module. The second sheet is the observed number of significant effects per
%   module. The third sheet is the expected number of non-significant effects
%   per module. The fourth sheet is the expected number of significant
%   effects per module.
%
%   module_mat_ratio -- a community X community X 2 matrix containing the
%   proportion of significant/non-significant observations to all observations per module. The first
%   sheet reflects the proportion of non-significant observations to all observations.
%   The second sheet represents the proportion of significant observations
%   to all observations.
%
%   module_mat_pvalue -- a community X community X 2 matrix containing the
%   pvalues and the thresholded significant effects via false discovery
%   rate. The first sheet contains the p-values. The second sheet is a
%   matrix of ones, zeros, and negative ones, which represent significantly greater than chance (1), not
%   significant (0), or significnatly less than chance (1).
%
%NOTES:
%  
%   The output matrices will be sorted in ascending numerical order. This
%   may not reflect the sorting you prefer, and you may have to re-sort the
%   matrix afterwards.
%
%SEE ALSO: CalculateChisquarePvalues

%VERSION HISTORY%
%
%VERSION 1.0 --4/30/14
%   Initialized and Documented by Eric Feczko
%
%
moduleval = unique(modules);
tempsize = size(moduleval);
if moduleval(1) == 0
    nummods = tempsize(1)-1;
else
    nummods = tempsize(1);
end
module_mat_count = zeros(nummods,nummods,4);
module_mat_chisquare = zeros(nummods,nummods);
module_mat_ratio = zeros(nummods,nummods,2);
mh = size(m,1);
mw = size(m,2);
for i = 1:nummods
    for j = 1:nummods
        vectorstep = 0;
        for k = 1:mh
            for l = 1:mw
                if moduleval(1) == 0
                    if modules(k) == moduleval(i+1) && modules(l) == moduleval(j+1)
                        if m(k,l) == 0
                            module_mat_count(i,j,1) = module_mat_count(i,j,1) + 1;
                        elseif m(k,l) == 1
                            module_mat_count(i,j,2) = module_mat_count(i,j,2) + 1;
                        end
                    end
                else
                    if modules(k) == moduleval(i) && modules(l) == moduleval(j)
                        if m(k,l) == 0
                            module_mat_count(i,j,1) = module_mat_count(i,j,1) + 1;
                        elseif m(k,l) == 1
                            module_mat_count(i,j,2) = module_mat_count(i,j,2) + 1;
                        end
                    end
                end
            end
        end
    end
end
module_mat_sum = module_mat_count(:,:,1) + module_mat_count(:,:,2);
module_mat_ratio(:,:,2) = module_mat_count(:,:,2)./module_mat_sum;
module_mat_ratio(:,:,1) = module_mat_count(:,:,1)./module_mat_sum;
module_mat_ratio_sig = module_mat_ratio(:,:,2);
module_mat_ratio_notsig = module_mat_ratio(:,:,1);
expected_ratio_sig = mean(mean(module_mat_ratio_sig));
expected_ratio_notsig = mean(mean(module_mat_ratio_notsig));
module_mat_count(:,:,3) = module_mat_sum .* expected_ratio_notsig;
module_mat_count(:,:,4) = module_mat_sum .* expected_ratio_sig;
for i = 1:nummods
    for j = 1:nummods
        observed_values = zeros(1,2);
        expected_values = zeros(1,2);
        observed_values(1,1) = module_mat_count(i,j,1);
        observed_values(1,2) = module_mat_count(i,j,2);
        expected_values(1,1) = module_mat_count(i,j,3);
        expected_values(1,2) = module_mat_count(i,j,4);
        module_mat_chisquare(i,j) = chi2(observed_values,expected_values);
    end
end
if nargin > 2
module_mat_pvalue = CalculateChisquarePvalues(module_mat_chisquare,module_mat_count,df,pvalr);
end
end