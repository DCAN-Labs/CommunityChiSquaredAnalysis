function chi2 = chi2(observed,expected)

% CHI2    Chi square statistic for contingency table.
%
% CHI2 = CHI2(OBSERVED)
% CHI2 = CHI2(OBSERVED,EXPECTED)
%
% Only the one-dimensional case supported so far
%
% See also CHI2PDF,CHI2INV, CROSSTAB, XTAB, BINOMIAL_LOGL.

% Original coding by Alex Petrov, Ohio State University
% $Revision: 1.0 $  $Date: 1999/12/18 17:20 $
%
% Part of the utils toolbox version 1.1 for MATLAB version 5 and up.
% http://alexpetrov.com/softw/utils/
% Copyright (c) Alexander Petrov 1999-2006, http://alexpetrov.com
% Please read the LICENSE and NO WARRANTY statement in ../utils_license.m

k = length(observed) ;     % assume one-dimensional vector
N = sum(observed) ;

if (nargin==1)
  expected = repmat(N/k,size(observed)) ;    % uniform default
end

chi2 = sum(((observed - expected) .^ 2) ./ expected) ;

%%%%% End of CHI2.M
