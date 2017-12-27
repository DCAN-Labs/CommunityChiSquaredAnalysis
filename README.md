##Community Chisquared Analysis (CCA) MATLAB Package##
###Created by Eric Feczko 4/30/14###

CCA comprises a set of matlab scripts designed to perform an exploratory
chi-squared test to identify interesting intra-/inter-network effects from a set
of mass univariate tests. The code generally requires two inputs determined
by the user:

1: A NxN matrix of binary values denoting significant effects; N refers to the
   number of nodes across all networks
2: A Nx1 matrix of numbers representing network assignments for each node.

Unfortunately, the code can only be run within a MATLAB session. Fortunately, 
no toolboxes should be required to run any commands within the package.

Below are the current commands I recommend calling when using the package:

[pvalues_modchisquare,observed_chisquare,observed_count,observed_ratio,...
...permuted_chisquare,permuted_count,permuted_ratio] = ...
...PermuteSignificantClusteringEffects(msig,modules,npermutations,pvalr): 
PermuteSignficantClusteringEffects will conduct the Chi-squared permutation test
on the input data. "msig" is the NxN matrix, "modules" is the Nx1 matrix.
The user can set the number of permutations "npermutations" to determine 
the p value. "pvalr" is the threshold used to determine the family-wise error
(FWE) rate via false-discovery rate (FDR) method 
(see: Benjamini and Hochberg, 1995).

PermuteCrossCorrelationMatrix(m,behavior,test_type,npermutations): 
PermuteCrossCorrelationMatrix is a slow but much more accurate method for
assessing the FWE for a Chi-squared analysis. Works for linear regression but 
not other statistical methods. "m" is a NxNxS data matrix where S is the number
of cases. "behavior" is a Nx1 matrix containing the outcome measure.

[observed_marginal_mean,observed_standard_error] = ...
...ProduceMarginalMeansFromChiSquared(msig,modules,pvalues_modchisquare,...
...varargin):
ProduceMarginalMeansFromChiSquared generates estimated marginal means for 
more complicated tests (e.g. ANOVAs). "pvalues_modchisquare" refers to the MxM
output matrix from PermuteSignificantClusteringEffects. "varargin" represents
multiple inputs to the function. Each group's data (a subset from "m") should be
entered in the same order as the "msig" and "modules" matrices.

SingleErrorPlot(observed_marginal_mean,observed_standard_error,filename,...
...colors,group_names):
SingleErrorPlot produces a plot of the marginal means for a single inter- or
intra-network comparison. Inputs come from the outputs of 
ProduceMarginalMeansFromChiSquared.

TrellisErrorPlot(observed_marginal_mean,observed_standard_error,...
...output_directory,colors):
TrellisErrorPlot produces a trellis plot of marginal means for all network 
comparisons. Inputs come from the outputs of ProduceMarginalMeansFromChiSquared.