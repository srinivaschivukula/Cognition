function Params=BasicClassifier

% Options For Collapsing Neural Data Through Time into a feature
Params.TimeCollapseMethod='Null';
Params.NumBins=nan; % Number of TimeSamples to Average across

% function to call to return Training Data
Params.TrainingDataFcn=[];

% Feature Reduction Options
% Params.FeatRedConfig=@Predictor.FeatRedConfig.BasicPropCutOff;
% Params.FeatRedConfig=[];
Params.FeatRedConfig=@Analyze.FaceScratch.MinimalReduction;
% Params.FeatRedConfig=[];

% Classifier Options
Params.ClassifierConfig.Type='LDA';
Params.ClassifierConfig.Arguements={'DiscrimType','diagLinear','FillCoeffs','off'};

Params.binSizeInSecs=0.05;