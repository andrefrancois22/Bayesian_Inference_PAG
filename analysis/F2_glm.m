clear all; close all; clc;
% ==> directories
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '../data/';

% ==> will store log likelihood ratios etc
llr = nan(1,29);

% ==> what is the session
for iSfln = 1:29 

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % params: matrix with the model parameters [trial x parameter]
    params = load([drc,'trial_DV_params_iS_',num2str(iSfln),'.mat']);
    params = params.ps_cat;
    fprintf('Finished loading matrix with the model parameters for iS = %d... \n',iSfln)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % dvs are model predicted DV trajectories [trial x time] (Category DVs)
    dvs = load([drc,'trial_DV_traj_iS_',num2str(iSfln),'.mat']);
    dvs = dvs.dv_cat;
    fprintf('Finished loading model predicted DV trajectories (Categorical DVs) for iS = %d... \n',iSfln)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % ==> load session data
    if iSfln > 9
        load([dataPath,'/dataSet_',num2str(iSfln),'.mat']);
    elseif iSfln <= 9
        load([dataPath,'/dataSet_0',num2str(iSfln),'.mat']);
    end

    % ==> choice
    cho = S.beh.choiceCat;
    % ==> context, contrast, orientation indicator variables
    ctx = S.exp.taskContext; ctr = S.exp.stimContrast; ori = S.exp.stimOriDeg;


    % ~~~~~~~~~~~~~~~~~~~~~~~~ Cat DV (signed) peak ~~~~~~~~~~~~~~~~~~~~~~~ %
    % ==> get max (maximum is either the positive peak if it's highest, 
    % or abs of negative peak if it's lowest)
    %mx = [max(dvs, [], 2), abs(min(dvs, [], 2))];
    %mxv = max(mx, [], 2);
    % => signed
    %mxv(mx(:,1) <= mx(:,2)) = -mxv(mx(:,1) <= mx(:,2));
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~ Cat DV dynamic range ~~~~~~~~~~~~~~~~~~~~~~~ %
    dynran = max([max(dvs, [], 2) - dvs(:,1), -(min(dvs, [], 2) - dvs(:,1))], [], 2);
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    
    % x axis is orientation (the values differ for FN and JP!). Use unique values
    or = unique(ori)'; cx = unique(ctx)'; cr = unique(ctr)';
    
    % ==> glm model design matrix
    X = [ori, dynran];
    y = cho;
    y(y == -1) = 0;
    
    % ==> logistic regression (ori and dynran)
    [b1,~,stats1] = glmfit(X,y,'binomial','link','logit');    
    % => yhat
    yh1 = glmval(b1,X,'logit');
    % ==> model (just orientation)
    [b0,~,stats0] = glmfit(X(:,1),y,'binomial','link','logit');    
    % => yhat
    yh0 = glmval(b0,X(:,1),'logit');
    
    % => acc (model 1)
    acc1 = sum(y==(yh1>=0.5))/length(y);
    % => acc (null model 0)
    acc0 = sum(y==(yh0>=0.5))/length(y);

    % ==> store model comparison metric
    llr(iSfln) = acc1/acc0;
    
%     figure(); plot(1:length(yh1(1:50)),yh1(1:50)>0.5,'r.-');
%     hold on;
%     plot(1:length(yh1(1:50)),y(1:50),'bo');
%     drawnow;
    
    
    % => console updates...
%     fprintf('Computed results for session %d of %d... \n',iSfln,29)
    fprintf('Model 1 Accuracy = %d...\n', acc1)
    fprintf('Model 0 Accuracy = %d...\n\n', acc0)
    
end

% ==> load dvCatperf data
dvCatPerf = load([drc,'dvCatPerf.mat'],'dvCatPerf');
dvCatPerf = dvCatPerf.dvCatPerf;
% => sort and return index
[~,idx] = sort(dvCatPerf);

figure();
plot(dvCatPerf(idx),llr(idx),'o');