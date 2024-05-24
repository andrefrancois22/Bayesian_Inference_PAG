clear all; close all; clc;
% ==> directories
dataPath     = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data';
drc = '../data/';
% ==> directory containing PF curve fit functions
pfc_functions_dr = '/home/thomas/Desktop/Bayesian_Inference_PAG/simulation/F1/pfc_functions/';
addpath(pfc_functions_dr)

%% ==> compute DV dynamic range medians for each stimulus type (29, 7 x 2 x 2) e.g. population x orientation x context by contrast


dynf = @(x) max([max(x, [], 2) - x(:,1), -(min(x, [], 2) - x(:,1))], [], 2);

% ==> DV dynamic range medians for each stimulus
DV_dynr_meds = nan(29,7,2,2);

% ==> offset Cat
ofs = [];
% ==> speed rise cat
sri = [];
% ==> dynr
dynr = [];

% ==> what is the session
for iS = 1:29 

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % params: matrix with the model parameters [trial x parameter]
    params = load([drc,'trial_DV_params_iS_',num2str(iS),'.mat']);
    params = params.ps_cat;
    fprintf('Finished loading matrix with the model parameters for iS = %d... \n',iS)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % dvs are model predicted DV trajectories [trial x time] (Category DVs)
    dvs = load([drc,'trial_DV_traj_iS_',num2str(iS),'.mat']);
    dvs = dvs.dv_cat;
    fprintf('Finished loading model predicted DV trajectories (Categorical DVs) for iS = %d... \n',iS)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % ==> store offset
    ofs = [ofs; mean(abs(params(:,1)))];
    
    % ==> store speed rise
    sri = [sri; mean(params(:,4))];
    
    % ==> dynamic range
    dynr = [dynr; mean(dynf(dvs))]; %params(2);
end

% ==> DV parameters
% offsetCat     = params(1);  % Bias in favor of CW or CCW choice
% scalarCat     = params(2);  % Controls dynamic range of cat dv
% spreadRiseCat = params(3);  % Controls speed rise of cat dv
% midRiseCat    = params(4);  % Controls half rise time of cat dv
% decayCstCat   = params(5);  % Controls decay after peak of cat dv
% offsetDir     = params(6);  % Bias in favor of LW or RW saccade
% scalarDir     = params(7);  % Controls dynamic range of dir dv
% spreadRiseDir = params(8);  % Controls speed rise of cat dv
% midRiseDir    = params(9) + midRiseCat - 1.65*spreadRiseCat + 1.65*spreadRiseDir; % Controls half rise time of dir dv

figure(1); set(gcf,'color','white');
subplot(1,2,1);
scatter(ofs,dynr,50,'ko','filled'); title(num2str(corr(ofs,dynr)));
xlabel('mean offset');
ylabel('mean dynamic range');
axis square;
subplot(1,2,2);
scatter(sri, dynr,50,'ko','filled'); title(num2str(corr(ofs,sri)));
xlabel('mean (half rise time)');
ylabel('mean dynamic range')
axis square;