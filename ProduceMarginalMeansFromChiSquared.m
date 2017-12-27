function [observed_marginal_mean,observed_standard_error] = ProduceMarginalMeansFromChiSquared(msig,modules,pvalues_modchisquare,varargin)
%ProduceMarginalMeansFromChiSquared will take the output from the
%chi-squared test and produce estimated marginal means and standard error
%for modules showing significant clusters
%
%USAGE: [observed_marginal_mean,observed_standard_error] = ProduceMarginalMeansFromChiSquared(datamat,modules,datamat_group1,datamat_group2,...);
%
%INPUTS:
%
%   msig -- A binarized matrix for signifiant (1) and not significant (0)
%   events. Each row/column in the matrix reflects a single region of interest (ROI). 
%   Typically, a binarized matrix results from thresholding a
%   set of significance tests (e.g. t-test or pearson's correlation)
%
%   modules -- A 1-D matrix where the length of the matrix is the same size
%   as the number of ROIs in the binarized matrix. Each item in the matrix
%   is a number that represents the community assignment for the given ROI.
%   The order of items should match the order of the binarized matrix.
%
%   pvalues_modchisquare -- The output from the chi-square analysis to
%   determine which modules need to be plotted. Only the second sheet is used here. 
%   The order for the pvalues_modchisquare should be numerical 
%   (from smallest to largest) and reflect the community assignments 
%   in modules.
%
%   datamat_groupx -- The raw data for each group should be provided as
%   separate variables at the end. The order for the correlation matrices
%   should be the same as the msig and modules order. 
%
%
%OUTPUTS:
%
%   Observed marginal means -- a community X community X N matrix containing the
%   marginal means for each group analyzed, where N is the total number of
%   groups. Marginal means are calculated by adjusting the grand mean
%
%   observed standard error -- a community X community X N matrix containig
%   the error for each group analyzed
%
%NOTES:
%  
%   The output matrices will be sorted in ascending numerical order. This
%   may not reflect the sorting you prefer, and you may have to re-sort the
%   matrix afterwards.
%
%SEE ALSO: PermuteSignificantClusteringEffects CountSignificantEffectsByModules

%VERSION HISTORY%
%
%VERSION 1.0 --2/7/17
%   -Initialized by Eric Feczko
%
%

%Declare variables and check for the number of communities and groups
ngroups = length(varargin);
all_data = varargin{1};
group_vector{1} = 1:size(varargin{1},3);
index = group_vector{1}(end)+1;
if ngroups > 1
    for i = 2:ngroups
        all_data(:,:,end+1:end+size(varargin{i},3)) = varargin{i};
        group_vector{i} = index:index+size(varargin{i},3)-1;
        index = group_vector{i}(end) + 1;
    end
end
ncommunities = length(unique(modules));
observed_marginal_mean = nan(ncommunities,ncommunities,ngroups);
observed_standard_error = observed_marginal_mean;
%cycle through communities to check for significant chi-squared findings,
for i = 1:ncommunities
    for j = i:ncommunities
        %check for significant clustering, if it exists start to collate
        %data from significant tests
        if pvalues_modchisquare(i,j) == 1
            for g=1:ngroups
                temp_all_data = all_data(modules == i,modules == j,group_vector{g});
                nsubs = size(temp_all_data,3);
                sigvals = repmat(msig(modules==i,modules==j),1,1,nsubs);
                temp_grand_mean = mean(temp_all_data(sigvals==1),'omitnan');
                index = 1;
                nsig = size(find(sigvals(:,:,1)==1),1);
                subject_adj_vector = zeros(nsig*size(temp_all_data,3),1);
                for k=1:size(temp_all_data,3)
                    subject_data = temp_all_data(:,:,k);
                    subject_vector = subject_data(sigvals(:,:,1)==1);
                    subject_mean = mean(subject_vector,'omitnan');
                    subject_adj_vector(index:index+length(subject_vector)-1,1) = subject_vector - subject_mean + temp_grand_mean;
                    index = size(subject_adj_vector,1)+1;
                end
                observed_marginal_mean(i,j,g) = mean(subject_adj_vector,'omitnan');
                observed_standard_error(i,j,g) = std(subject_adj_vector,'omitnan')/sqrt(size(subject_adj_vector,1));
                observed_marginal_mean(j,i,g) = observed_marginal_mean(i,j,g);
                observed_standard_error(j,i,g) = observed_standard_error(i,j,g);
                clear subject_adj_vector
            end
        end
    end
end