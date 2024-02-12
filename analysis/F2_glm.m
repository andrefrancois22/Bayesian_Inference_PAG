clear all; close all; clc;
% ==> directories
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '../data/';

% ==> for storing 4 model design matrices
% => context 1 & contrast 1
cx1_cr1_X = cell([29,1]);
% => context 1 & contrast 2
cx1_cr2_X = cell([29,1]);

% => context 2 & contrast 1
cx2_cr1_X = cell([29,1]);
% => context 2 & contrast 2
cx2_cr2_X = cell([29,1]);



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
    or = unique(ori)'; 
    cx = unique(ctx)'; 
    cr = unique(ctr)';
    

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % => index context 1, and contrast 1 (low contrast)
    idx11 = (ctx == 1) & (ctr == cr(1));
    % 2) ==> dummy coding orientation    
       
    % => context 1 & contrast 1
    cx1_cr1_X{iSfln} = [ori(idx11), dynran(idx11)];
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % => index context 1, and contrast 2 (hi contrast)
    idx12 = (ctx == 1) & (ctr == cr(2));
    % 2) ==> dummy coding orientation        
       
    % => context 1 & contrast 1
    cx1_cr2_X{iSfln} = [ori(idx12), dynran(idx12)];    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % => index context 2, and contrast 1 (low contrast)
    idx21 = (ctx == -1) & (ctr == cr(1));
    % 2) ==> dummy coding orientation        
       
    % => context 1 & contrast 1
    cx2_cr1_X{iSfln} = [ori(idx21), dynran(idx21)];
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % => index context 2, and contrast 2 (low contrast)
    idx22 = (ctx == -1) & (ctr == cr(2));
    % 2) ==> dummy coding orientation        
       
    % => context 1 & contrast 1
    cx2_cr2_X{iSfln} = [ori(idx22), dynran(idx22)];    
    % ==> behavioral choice outcome is y
    y22 = cho(idx22);
    y22(y22 == -1) = 0;
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
    
end
 
% % ==> logistic regression (ori and dynran)
% [b1,~,stats1] = glmfit(X,y,'binomial','link','logit');    
% % => yhat
% yh1 = glmval(b1,X,'logit');
% % ==> model (just orientation)
% [b0,~,stats0] = glmfit(X(:,1),y,'binomial','link','logit');    
% % => yhat
% yh0 = glmval(b0,X(:,1),'logit');

% ==> load dvCatperf data
dvCatPerf = load([drc,'dvCatPerf.mat'],'dvCatPerf');
dvCatPerf = dvCatPerf.dvCatPerf;
% => sort and return index
[~,idx] = sort(dvCatPerf);

figure();
plot(dvCatPerf(idx),llr(idx),'o');



% figure(); plot(1:length(yh1(1:50)),yh1(1:50)>0.5,'r.-');
% hold on;
% plot(1:length(yh1(1:50)),y(1:50),'bo');
% drawnow;