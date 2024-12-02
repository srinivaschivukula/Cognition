function Out=dpcaAnalysis(firingRatesAverage,firingRates,trialNum,time,varargin)

% defaults
dPCA_noReg=1;
dPCA_Reg=1;
% lambdas=1e-07 * 1.5.^[0:1:20];% allows less regularized
% lambdas=1e-07 * 1.5.^[0:1:30]; % allows even more regularized
lambdas=1e-07 * 1.5.^[0:1:25];% default
% lambdas=1e-07 * 1.5.^[0:1:23];% default
% Time events of interest (e.g. stimulus onset/offset, cues etc.)
[varargin,timeEvents]=Utilities.ProcVarargin(varargin,'timeEvents',0);
[varargin,plotPCAMarg]=Utilities.ProcVarargin(varargin,'plotPCAMarg',1);
[varargin,plotPCABasic]=Utilities.ProcVarargin(varargin,'plotPCABasic',1);
% how many cross-validation iterations to perform in estimating lambda
[varargin,numRepReg]=Utilities.ProcVarargin(varargin,'numRepReg',10);

% number of dpca components
[varargin,numComponents]=Utilities.ProcVarargin(varargin,'numComponents',20);
% number of shuffles for significance testing
[varargin,numShuffles]=Utilities.ProcVarargin(varargin,'numShuffles',100);
% how many cross-validation iterations to perform for decode
[varargin,numRepDecode]=Utilities.ProcVarargin(varargin,'numRepDecode',100);

[varargin,dPCA_noReg]=Utilities.ProcVarargin(varargin,'dPCA_noReg',dPCA_noReg);
[varargin,dPCA_Reg]=Utilities.ProcVarargin(varargin,'dPCA_Reg',dPCA_Reg);

[varargin,decode_noReg]=Utilities.ProcVarargin(varargin,'decode_noReg');
[varargin,decode_Reg]=Utilities.ProcVarargin(varargin,'decode_Reg');

% name to associate with each
[varargin,namePrefix]=Utilities.ProcVarargin(varargin,'namePrefix','');
[varargin,baseDirectory]=Utilities.ProcVarargin(varargin,'baseDirectory','');
figsDirectory=fullfile(baseDirectory,'Figs'); if ~exist(figsDirectory); mkdir(figsDirectory); end
dataDirectory=fullfile(baseDirectory,'Data');if ~exist(dataDirectory); mkdir(dataDirectory); end

Utilities.ProcVarargin(varargin);

N = size(firingRatesAverage,1);    % number of neurons
T = size(firingRatesAverage,ndims(firingRatesAverage));     % number of time points
S = size(firingRatesAverage,2);      % body-parts (cheek, shoulder)
D = size(firingRatesAverage,3);  % person (nancy,tyson)
E = size(firingRates,ndims(firingRates));     % maximal number of trial repetitions

% trialNum: N x S x D
% firingRates: N x S x D x T x maxTrialNum
% firingRatesAverage: N x S x D x T
%
% N is the number of neurons
% S is the number of stimuli conditions (F1 frequencies in Romo's task)
% D is the number of decisions (D=2)
% T is the number of time-points (note that all the trials should have the
% same length in time!)
%
% trialNum -- number of trials for each neuron in each S,D condition (is
% usually different for different conditions and different sessions)
%
% firingRates -- all single-trial data together, massive array. Here
% maxTrialNum is the maximum value in trialNum. E.g. if the number of
% trials per condition varied between 1 and 20, then maxTrialNum = 20. For
% the neurons and conditions with less trials, fill remaining entries in
% firingRates with zeros or nans.
%
% firingRatesAverage -- average of firingRates over trials (5th dimension).
% If the firingRates is filled up with nans, then it's simply
%    firingRatesAverage = nanmean(firingRates,5)
% If it's filled up with zeros (as is convenient if it's stored on hard
% drive as a sparse matrix), then
%    firingRatesAverage = bsxfun(@times, mean(firingRates,5), size(firingRates,5)./trialNum)

%% Define parameter grouping

% parameter groupings
% 1 - body-part
% 2 - person
% 3 - time
% There are three pairwise interactions:
%    [1 3] - bodypart/time interaction
%    [2 3] - person/time interaction
%    [1 2] - bodypart/person interaction
% And one three-way interaction:
%    [1 2 3] - rest

% Here we group stimulus with stimulus/time interaction etc. Don't change
% that if you don't know what you are doing

% combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}};
combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
timeMarginalization=2;
margNames = {'Body-Part','Person','Time','Interaction'};
% margNames = {'Body-Part','Person','Time'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% combinedParams =  {{1, [1 2]}, {2}};
% timeMarginalization=2;
% margNames = {'Stimulus','Condition-Independent'};
% margColours = [23 100 171; 187 20 25]/256;

% % setting random number of repetitions for each neuron and condition
ifSimultaneousRecording = true;  % change this to simulate simultaneous 
                                 % recordings (they imply the same number 
                                 % of trials for each neuron)



%% Step 1: PCA of the dataset
if plotPCABasic
    X = firingRatesAverage(:,:);
    X = bsxfun(@minus, X, mean(X,2));
    [W,~,~] = svd(X);
    
    
    % computing explained variance
    explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
        'combinedParams', combinedParams);
    
    % a bit more informative plotting
    dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
        'explainedVar', explVar, ...
        'marginalizationNames', margNames, ...
        'time', time,                        ...
            'timeEvents', timeEvents,               ...
        'marginalizationColours', margColours);
    
    Analyze.AddFigLabel('',[],'PCA');
end
%% Step 2: PCA in each marginalization separately
if plotPCAMarg
    dpca_perMarginalization(firingRatesAverage, @dpca_plot_default, ...
        'combinedParams', combinedParams);
    Analyze.AddFigLabel('',[],'PCA-Marginalized');
end

%% Step 3: dPCA without regularization

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization
if dPCA_noReg
    [W,V,whichMarg] = dpca(firingRatesAverage, numComponents, ...
        'combinedParams', combinedParams);
    
    explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
        'combinedParams', combinedParams);
    Out.noReg.W=W;
    Out.noReg.V=V;
    Out.noReg.explVar=explVar;
    Out.noReg.margNames=margNames;
    Out.noReg.whichMarg=whichMarg;
    
    dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
        'explainedVar', explVar, ...
        'marginalizationNames', margNames, ...
        'marginalizationColours', margColours, ...
        'whichMarg', whichMarg,                 ...
        'time', time,                        ...
        'timeEvents', timeEvents,               ...
        'timeMarginalization', timeMarginalization, ...
        'legendSubplot', 16);
    
    Analyze.AddFigLabel('',[],'dPCA-NoReg');
end

%% Step 4: dPCA with regularization

% This function takes some minutes to run. It will save the computations
% in a .mat file with a given name. Once computed, you can simply load
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')



if dPCA_Reg
    optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
        'combinedParams', combinedParams, ...
        'numComps', numComponents,...
        'lambdas', lambdas,...
        'numRep', 10);  % increase this number to ~10 for better accuracy
        
    Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

    [W,V,whichMarg] = dpca(firingRatesAverage, numComponents, ...
        'combinedParams', combinedParams, ...
        'lambda', optimalLambda, ...
        'Cnoise', Cnoise);
    
    explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
        'combinedParams', combinedParams);
    
    Out.Reg.W=W;
    Out.Reg.V=V;
    Out.Reg.explVar=explVar;
    Out.Reg.margNames=margNames;
    Out.Reg.whichMarg=whichMarg;
    
    dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
        'explainedVar', explVar, ...
        'marginalizationNames', margNames, ...
        'marginalizationColours', margColours, ...
        'whichMarg', whichMarg,                 ...
        'time', time,                        ...
        'timeEvents', timeEvents,               ...
        'timeMarginalization', timeMarginalization,           ...
        'legendSubplot', 16);
    
    Analyze.AddFigLabel('',[],'dPCA-Regularized');
end

%% Optional: estimating "signal variance"

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams, ...
    'Cnoise', Cnoise, 'numOfTrials', trialNum);

% Note how the pie chart changes relative to the previous figure.
% That is because it is displaying percentages of (estimated) signal PSTH
% variances, not total PSTH variances. See paper for more details.

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', timeMarginalization,           ...
    'legendSubplot', 16);

%% Decoding

if decode_Reg && dPCA_Reg
    
    decodingClasses = {[(1:S)' (1:S)'], repmat([1:2], [S 1])};
    accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
        'lambda', optimalLambda, ...
        'combinedParams', combinedParams, ...
        'decodingClasses', decodingClasses, ...
        'numRep', numRepDecode, ...        % increase to 100
        'filename', 'tmp_classification_accuracy.mat');
    %%
    accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
        'lambda', optimalLambda, ...
        'combinedParams', combinedParams, ...
        'decodingClasses', decodingClasses, ...
        'numRep', numRepDecode, ...        % increase to 100
        'numShuffles', numShuffles, ...  % increase to 100 (takes a lot of time)
        'filename', 'tmp_classification_accuracy.mat');
    
    componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
    
    dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
        'explainedVar', explVar, ...
        'marginalizationNames', margNames, ...
        'marginalizationColours', margColours, ...
        'whichMarg', whichMarg,                 ...
        'time', time,                        ...
        'timeEvents', timeEvents,               ...
        'timeMarginalization', timeMarginalization,           ...
        'legendSubplot', 16,                ...
        'componentsSignif', componentsSignif);
end
