function Params=MinimalReduction

Params.Tests.usePCA=0;
Params.nPrinComp=20;
Params.plotPCA=0;
Params.robustPCA=0;

% Per feature stat test.
Params.Tests.anovaTest=0;
Params.anovaThresh=.05;

% Data Normalization
Params.Tests.zscore=0;


% Rate Restrictions
Params.Tests.DeflectionCoeffTest=0;
% Params.minFiringRate=.25;

% Rate Restrictions
Params.Tests.FiringRateTest=0;
Params.minFiringRate=.25;

% at least two categories must be seperate by a firing rate greater
% than
Params.Tests.minCategoryRateDifferenceTest=0;
Params.maxCatFRDiff_thresh=1;

% at least two categories must be seperate by a dPrime greater
% than
Params.Tests.dPrimeTest=0;
Params.maxCatDPrime_thresh=.25;

% At least one category must have a mean rate greater than
Params.Tests.minCategoryRateTest=0;
Params.maxMeanRate_thresh=.25;

% At least one category must have a median rate greater than
Params.Tests.minCategoryMedianTest=0;
Params.maxMedianRate_thresh=0;

% At least one category must have a median rate greater than
Params.Tests.useQualityThresh=0;
Params.QualityThreshold=4;

% RRelieff Algorithm
Params.Tests.useRelief=0; 