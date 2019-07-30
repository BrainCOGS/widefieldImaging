function  [coeff, FitInfo] = ridgeRegression(X, y, lambdaRatio, numCVFolds, numCVSamples, jitter) 

% [coeff, fitInfo] = ridgeRegression(X, y, lambdaRatop, numCVFolds, numCVSamples) 
% wrapper to perform ridge regression using a CV object and having similar
% output structure to lasso. Calls the matlab function ridge.m
% X is predictor matrix, y is response vector, lambdaRatio is vector
% containing the lambda hyperparameter. numCVfolds and numCVsamples are to
% set up k-fold CV object. alternatively numCVFolds can be the CV object
% itself
% coeff is n params x n lambda matrix for full data solution, fitInfo
% contains CV info similar to lasso and elasticNetRegression

 %% Default arguments
if nargin < 3 || isempty(lambdaRatio)
  lambdaRatio           = flip(10.^(linspace(-3, log10(1), 30)'));
else
  lambdaRatio           = sort(lambdaRatio(:), 'descend');
end
if nargin < 4 || isempty(numCVFolds)
  numCVFolds            = 3;
end
if nargin < 5 || isempty(numCVSamples)
  numCVSamples          = 10;
end
if nargin < 6 || isempty(jitter)
  jitter                = [1e-6 1e-5];
end

%% Cross-validation setup
if isnumeric(numCVFolds)
  warning('off', 'stats:cvpartition:KFoldMissingGrp');
  pseudoExp           = cvpartition(y, 'KFold', numCVFolds);
  for iMC = 2:numCVSamples
    pseudoExp(iMC)    = repartition(pseudoExp(1));
  end
  warning('on', 'stats:cvpartition:KFoldMissingGrp');
else
  pseudoExp           = numCVFolds;
  numCVFolds          = pseudoExp(1).NumTestSets;
end

%% central value fits
FitInfo.CVExperiments     = pseudoExp;
FitInfo.Coefficients      = ridge(y, X, lambdaRatio);
coeff                     = FitInfo.Coefficients;

%% Cross validated goodness-of-fit
mse                       = nan(numel(lambdaRatio), numel(pseudoExp), numCVFolds);
rSquared                  = nan(numel(lambdaRatio), numel(pseudoExp), numCVFolds);
FitInfo.CVTrainPrediction = cell(numel(pseudoExp), numCVFolds);
FitInfo.CVTestPrediction  = cell(numel(pseudoExp), numCVFolds);

for iMC = 1:numel(pseudoExp)
  for iFold = 1:numCVFolds
    % Select train/test sets with user-specified pruning if so desired
    if isfield(pseudoExp(iMC).training(iFold),'idx')
      iTrain              = pseudoExp(iMC).training(iFold).idx;
      iTest               = pseudoExp(iMC).test(iFold).idx;
    else
      iTrain              = pseudoExp(iMC).training(iFold);
      iTest               = pseudoExp(iMC).test(iFold);
    end
    % Perform fit for this CV experiment and generate/test prediction
    mcCoeff               = ridge(y(iTrain), X(iTrain,:), lambdaRatio, 0); % zero also returns bias term

    % train prediction
    FitInfo.CVTrainPrediction{iMC,iFold} = ridgePrediction(X(iTrain,:),mcCoeff);
    [FitInfo.CVTestPrediction{iMC,iFold}, mse(:,iMC,iFold), rSquared(:,iMC,iFold)] ...
                          = ridgePrediction(X(iTest,:), mcCoeff, y(iTest), jitter);
  end
end

%% Collect output in same format as Matlab's lasso()
FitInfo.Lambda            = lambdaRatio;
FitInfo.DF                = sum(coeff ~= 0, 1);
FitInfo.MSE               = mean(mse(:,:), 2, 'omitnan')';
FitInfo.SE                = std(mse(:,:), 0, 2, 'omitnan')';
FitInfo.RSquared          = mean(rSquared(:,:), 2, 'omitnan')';
FitInfo.RSquaredErr       = std(rSquared(:,:), 0, 2, 'omitnan')';
[~, FitInfo.IndexMinMSE]  = min(FitInfo.MSE);
FitInfo.LambdaMinMSE      = FitInfo.Lambda(FitInfo.IndexMinMSE);
FitInfo.Index1SE          = find(   FitInfo.MSE(1:FitInfo.IndexMinMSE)    ...
                                <=  FitInfo.MSE(FitInfo.IndexMinMSE)      ...
                                  + FitInfo.SE(FitInfo.IndexMinMSE)       ...
                                , 1, 'first'                              ...
                                );
if isempty(FitInfo.Index1SE)
  FitInfo.Index1SE        = FitInfo.IndexMinMSE;
  FitInfo.Lambda1SE       = FitInfo.LambdaMinMSE;
else
  FitInfo.Lambda1SE       = FitInfo.Lambda(FitInfo.Index1SE);
end

end

%% predict data for different lambdas
function [yhat, mse, rSquared] = ridgePrediction(X, coeff, y, jitter)

if nargin < 3; y      = []; end
if nargin < 4; jitter = []; end

% add bias column
if size(X,2) < size(coeff,1)
  X = [ones(size(X,1),1) X];
end

nlambda = size(coeff,2);
npoints = size(X,1);
yhat    = zeros(npoints,nlambda);
for iL = 1:nlambda
  yhat(:,iL) = X*coeff(:,iL);
end

if isempty(y)
  mse      = [];
  rSquared = [];
else
  y        = repmat(y,[1 nlambda]);
  sse      = bsxfun(@minus, y, yhat).^2;
  sse      = sum(sse, 1);
  mse      = (sse ./ npoints)';
  muY      = repmat(mean(y),[npoints 1]);
  ssy      = sum((y - muY).^2);
  rSquared = (1 - sse ./ ssy)';
end

if ~isempty(jitter)
  yhat     = yhat + jitter(1)*randn(size(yhat));% + rand(size(yhat)) * (jitter(2) - jitter(1));
end

end
