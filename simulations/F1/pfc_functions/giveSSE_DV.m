%% giveSSE_DV
% This script returns the SSE of a descriptive model of the dynamic 
% evolution of a population decision variable. Depending on the parameter
% bounds, this model can capture different neural coding strategies : the
% abstract hypothesis, the intentional hypothesis, and the mixture
% hypothesis.

function [SSE, dvCatPred, dvDirPred, dvCatPredHR, dvDirPredHR, timeSacHR] = giveSSE_DV(params, timeSac, dvCat, dvDir)

% Decode function arguments
offsetCat     = params(1);  % Bias in favor of CW or CCW choice
scalarCat     = params(2);  % Controls dynamic range of cat dv
spreadRiseCat = params(3);  % Controls speed rise of cat dv
midRiseCat    = params(4);  % Controls half rise time of cat dv
decayCstCat   = params(5);  % Controls decay after peak of cat dv
offsetDir     = params(6);  % Bias in favor of LW or RW saccade
scalarDir     = params(7);  % Controls dynamic range of dir dv
spreadRiseDir = params(8);  % Controls speed rise of cat dv
midRiseDir    = params(9) + midRiseCat - 1.65*spreadRiseCat + 1.65*spreadRiseDir; % Controls half rise time of dir dv

% Compute model predictions (trajectory is shaped like cumulative Gaussian...)
dvCatPred = offsetCat + scalarCat * normcdf(timeSac', midRiseCat, spreadRiseCat);
dvDirPred = offsetDir + scalarDir * normcdf(timeSac', midRiseDir, spreadRiseDir);

% ...Followed by decay in cat dv
decayTime         = midRiseCat + 2.5*spreadRiseCat;
selInd            = timeSac' > decayTime;
dvCatPred(selInd) = dvCatPred(selInd) .* ((1 - decayCstCat) .^ (.001*(timeSac(selInd)' - decayTime)));

% Compute high resolution model prediction
timeSacHR             = linspace(timeSac(1), timeSac(end), 200)';
dvCatPredHR           = offsetCat + scalarCat * normcdf(timeSacHR, midRiseCat, spreadRiseCat);
dvDirPredHR           = offsetDir + scalarDir * normcdf(timeSacHR, midRiseDir, spreadRiseDir);
selIndHR              = timeSacHR > decayTime;
dvCatPredHR(selIndHR) = dvCatPredHR(selIndHR) .* ((1 - decayCstCat) .^ (.001*(timeSacHR(selIndHR) - decayTime)));

% Compute SSE
SSE = sum((dvCat - dvCatPred).^2 + (dvDir - dvDirPred).^2);
end

%%
