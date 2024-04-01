clear all; close all; clc;
% ==> directories
dataPath     = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data';
drc = '../data/';

% ==> directory containing PF curve fit functions
pfc_functions_dr = '/home/thomas/Desktop/Bayesian_Inference_PAG/simulation/pfc_functions/';
addpath(pfc_functions_dr)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ==> PF curve estimation Set options
fitFlag         = 1;       % 0 or 1, where 1 means do the fit right now
nBootStraps     = 100;     % determines number of non parametric bootstraps if fitFlag == 1
options         = optimoptions('fmincon');
options.Display = 'off';
% Set bounds on model parameters
% Model fit to choice data split by task context (1 & 2: guess rate; 3: perceptual uncertainty; 4 & 5: decision criterion)
startVec_M1 = [0.04 0.04 1 0 0];
% startVec_V2 = [0.04 0.04 3 -15 15];
LB_M1(1,1)  = 0;                         UB_M1(1,1) = 0.05;      % lapse rate
LB_M1(2,1)  = 0;                         UB_M1(2,1) = 0.05;      % lapse rate
LB_M1(3,1)  = 0.1;                       UB_M1(3,1) = 10;         % perceptual uncertainty
LB_M1(4,1)  = -200;                      UB_M1(4,1) = 200;        % decision criterion (-10->10 range originally)
LB_M1(5,1)  = -200;                      UB_M1(5,1) = 200;        % decision criterion
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%%
close all; clc;

% ==> dynamic range and other data
%res = {};

dbias = nan(29,1);
dpu   = nan(29,1);

figure(4); set(gcf,'color','white');

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
    timeSac = S.dec.timeSac;
    % Set time boundaries
    fitTimeSacBegin = -800;  % in ms, relative to saccade onset
    fitTimeSacEnd   = -50;
    [~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
    [~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));

    % time intervals
    t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));

    % ==> choice
    cho = S.beh.choiceCat;
    % ==> context, contrast, orientation indicator variables
    ctx = S.exp.taskContext; ctr = S.exp.stimContrast; ori = S.exp.stimOriDeg;
    % x axis is orientation (the values differ for FN and JP!). Use unique values
    or = unique(ori)'; cx = unique(ctx)'; cr = unique(ctr)';
    
    % ==> dynamic range for session
    dynr = max([max(dvs, [], 2) - dvs(:,1), -(min(dvs, [], 2) - dvs(:,1))], [], 2);
    
    % ==> dynamic range median split
    % ==> hi & lo dynr and hi & lo contrast and cw vs. ccw context
    
    % ==> store, for each session:
    % ---> dynamic range
    % ---> context indicator variable
    % ---> contrast indicator variable
    % ---> orientation
    % ---> choice
    %res{iSfln} = [dynr,ctx,ctr,ori,cho];
       
    % => indexing variables - dynamic range
    lo_dynr = (dynr(:,1) <  median(dynr(:,1)));    
    hi_dynr = (dynr(:,1) >= median(dynr(:,1)));

    % ==> cw context (congruent choices)
    % ==> cw context & lo contrast, and lo dynamic range
    cw_lcr_ldr_cng = (ctx == 1) & (ctr == min(ctr)) & (lo_dynr) & (cho == 1);    
    % ==> cw context & lo contrast, and hi dynamic range
    cw_lcr_hdr_cng = (ctx == 1) & (ctr == min(ctr)) & (hi_dynr) & (cho == 1);    
    % ==> ccw context (congruent choices)
    % ==> ccw context & lo contrast, and lo dynamic range
    ccw_lcr_ldr_cng = (ctx == -1) & (ctr == min(ctr)) & (lo_dynr) & (cho == -1);    
    % ==> ccw context & lo contrast, and hi dynamic range
    ccw_lcr_hdr_cng = (ctx == -1) & (ctr == min(ctr)) & (hi_dynr) & (cho == -1);      
        
    
    % ==> cw context (incongruent choices)
    % ==> cw context & lo contrast, and lo dynamic range
    cw_lcr_ldr_icg = (ctx == 1) & (ctr == min(ctr)) & (lo_dynr) & (cho == -1);    
    % ==> cw context & lo contrast, and hi dynamic range
    cw_lcr_hdr_icg = (ctx == 1) & (ctr == min(ctr)) & (hi_dynr) & (cho == -1);            
    % ==> ccw context (incongruent choices)
    % ==> ccw context & lo contrast, and lo dynamic range
    ccw_lcr_ldr_icg = (ctx == -1) & (ctr == min(ctr)) & (lo_dynr) & (cho == 1);    
    % ==> ccw context & lo contrast, and hi dynamic range
    ccw_lcr_hdr_icg = (ctx == -1) & (ctr == min(ctr)) & (hi_dynr) & (cho == 1);       
    
    % ==> ratio
    cho_cw_ct_cw_cr_lo_dynr_lo   = nan(1,7);
    cho_cw_ct_ccw_cr_lo_dynr_lo  = nan(1,7);
    
    cho_ccw_ct_cw_cr_lo_dynr_lo  = nan(1,7);
    cho_ccw_ct_ccw_cr_lo_dynr_lo = nan(1,7);
    
    cho_cw_ct_cw_cr_lo_dynr_hi   = nan(1,7);
    cho_cw_ct_ccw_cr_lo_dynr_hi  = nan(1,7);
    
    cho_ccw_ct_cw_cr_lo_dynr_hi  = nan(1,7);
    cho_ccw_ct_ccw_cr_lo_dynr_hi = nan(1,7);
    
    % ==> compute counts for each stimulus orientation
    for theta = 1:7
        % ==> low dynamic range split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % ==> cw congruent choice, lo ctr, lo dynr
        cho_cw_ct_cw_cr_lo_dynr_lo(theta)   = sum(cw_lcr_ldr_cng   & ori==or(theta));% / sum((ctx == 1) & (ctr == min(ctr)) & (lo_dynr)  & ori==or(theta));
        % ==> cw incongruent choice, lo ctr, lo dynr
        cho_cw_ct_ccw_cr_lo_dynr_lo(theta)  = sum(ccw_lcr_ldr_icg  & ori==or(theta));% / sum((ctx == -1) & (ctr == min(ctr)) & (lo_dynr) & ori==or(theta));
        
        % ==> ccw incongruent choice, lo ctr, lo dynr
        cho_ccw_ct_cw_cr_lo_dynr_lo(theta)  = sum(cw_lcr_ldr_icg   & ori==or(theta));% / sum((ctx == 1) & (ctr == min(ctr)) & (lo_dynr) & ori==or(theta));
        % ==> ccw congruent choice, lo ctr, lo dynr
        cho_ccw_ct_ccw_cr_lo_dynr_lo(theta) = sum(ccw_lcr_ldr_cng  & ori==or(theta));% / sum((ctx == -1) & (ctr == min(ctr)) & (lo_dynr) & ori==or(theta));
        
 
        % ==> high dynamic range split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         % ==> cw congruent choice, lo ctr, lo dynr
        cho_cw_ct_cw_cr_lo_dynr_hi(theta)   = sum(cw_lcr_hdr_cng   & ori==or(theta));% / sum((ctx == 1) & (ctr == min(ctr)) & (lo_dynr)  & ori==or(theta));
        % ==> cw incongruent choice, lo ctr, lo dynr
        cho_cw_ct_ccw_cr_lo_dynr_hi(theta)  = sum(ccw_lcr_hdr_icg  & ori==or(theta));% / sum((ctx == -1) & (ctr == min(ctr)) & (lo_dynr) & ori==or(theta));
        
        % ==> ccw incongruent choice, lo ctr, lo dynr
        cho_ccw_ct_cw_cr_lo_dynr_hi(theta)  = sum(cw_lcr_hdr_icg   & ori==or(theta));% / sum((ctx == 1) & (ctr == min(ctr)) & (lo_dynr) & ori==or(theta));
        % ==> ccw congruent choice, lo ctr, lo dynr
        cho_ccw_ct_ccw_cr_lo_dynr_hi(theta) = sum(ccw_lcr_hdr_cng  & ori==or(theta));% / sum((ctx == -1) & (ctr == min(ctr)) & (lo_dynr) & ori==or(theta));        
    end
    
%     % ~~~~~~~~ lo contrast & lo dynamic range ~~~~~~~~~
%     %   ==> cw  context ccw response (incongruent) 
%     % & ==> ccw context ccw response (congruent)
%     cho_ccw_ct_cw_cr_lo_dynr_lo  % = cw_lcr_ldr_icg;
%     cho_ccw_ct_ccw_cr_lo_dynr_lo % = ccw_lcr_ldr_cng;
% 
%     %   ==> cw  context cw response (congruent)    
%     % & ==> ccw context cw response (incongruent)
%     cho_cw_ct_cw_cr_lo_dynr_lo  % = cw_lcr_ldr_cng;
%     cho_cw_ct_ccw_cr_lo_dynr_lo % = ccw_lcr_ldr_icg;
%     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
%     % ~~~~~~~~ lo contrast & hi dynamic range ~~~~~~~~~
%     %   ==> cw  context ccw response (incongruent) 
%     % & ==> ccw context ccw response (congruent)
%     cho_ccw_ct_cw_cr_lo_dynr_hi  % = cw_lcr_hdr_icg;
%     cho_ccw_ct_ccw_cr_lo_dynr_hi % = ccw_lcr_hdr_cng;
% 
%     %   ==> cw  context cw response (congruent)    
%     % & ==> ccw context cw response (incongruent)
%     cho_cw_ct_cw_cr_lo_dynr_hi  % = cw_lcr_hdr_cng;
%     cho_cw_ct_ccw_cr_lo_dynr_hi % = ccw_lcr_hdr_icg;
%     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    
    % ==> PF curve fit
    fprintf('computing psychometric functions...\n')

    % ==> low dynamic range cases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ==> counts
    % ==> context 1 ccw response (incongruent) & ==> context 2 ccw response (congruent)
    nccw_dynr_lo = [cho_ccw_ct_cw_cr_lo_dynr_lo; cho_ccw_ct_ccw_cr_lo_dynr_lo]; 
    % ==> context 1 cw response (congruent)    & ==> context 2 cw response (incongruent)
    ncw_dynr_lo  = [cho_cw_ct_cw_cr_lo_dynr_lo;  cho_cw_ct_ccw_cr_lo_dynr_lo]; 
    % ==> pass through PF function fit
    obFun    = @(paramVec) giveNLL(paramVec, or, nccw_dynr_lo, ncw_dynr_lo, 'M1');    
    % ==> run 
    paramEst_M1_dynr_lo = fmincon(obFun, startVec_M1, [], [], [], [], LB_M1, UB_M1, [], options);
    % ==> return pf
    [~, predPF_dynr_lo] = giveNLL(paramEst_M1_dynr_lo, or, nccw_dynr_lo, ncw_dynr_lo, 'M1');
    %  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    % ==> high dynamic range cases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ==> counts 
    % ==> context 1 ccw response (incongruent) & ==> context 2 ccw response (congruent)
    nccw_dynr_hi = [cho_ccw_ct_cw_cr_lo_dynr_hi; cho_ccw_ct_ccw_cr_lo_dynr_hi]; 
    % ==> context 1 cw response (congruent)    & ==> context 2 cw response (incongruent)
    ncw_dynr_hi  = [cho_cw_ct_cw_cr_lo_dynr_hi;  cho_cw_ct_ccw_cr_lo_dynr_hi]; 
    % ==> pass through PF function fit
    obFun    = @(paramVec) giveNLL(paramVec, or, nccw_dynr_hi, ncw_dynr_hi, 'M1');    
    % ==> run 
    paramEst_M1_dynr_hi = fmincon(obFun, startVec_M1, [], [], [], [], LB_M1, UB_M1, [], options);
    % ==> return pf
    [~, predPF_dynr_hi] = giveNLL(paramEst_M1_dynr_hi, or, nccw_dynr_hi, ncw_dynr_hi, 'M1');
    %  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    % ==> plots
    subplot(5,6,iSfln); title(['Session: ', num2str(iSfln)]);
    plot(or,predPF_dynr_lo(1,:),'b.-')
    hold on; hold all;
    plot(or,predPF_dynr_lo(2,:),'r.-')
    hold on; hold all;
    plot(or,predPF_dynr_hi(1,:),'b+:')
    hold on; hold all;
    plot(or,predPF_dynr_hi(2,:),'r+:')    
    xlabel('stimulus orientation')
    ylabel('Proportion cw PF fit')
    drawnow;
    
    % ==> delta bias (this is a different computation from original)
    dbias(iSfln) = paramEst_M1_dynr_lo(5) - paramEst_M1_dynr_hi(5);    
    % ==> delta perceptual uncertainty
    dpu(iSfln)   = paramEst_M1_dynr_lo(3) - paramEst_M1_dynr_hi(3);
end

% ==> pearson correlation between delta bias and delta uncertainty
% ==> ommitting session 1 here!
[r,p] = corr(dbias(2:end),dpu(2:end));

figure(); set(gcf,'color','white');
plot(dbias(2:end),dpu(2:end),'k+')
title(['r = ',num2str(r), ', p = ',num2str(p)]);
xlabel('\Delta bias')
ylabel('\Delta perceptual uncertainty (slope)')

%% ==> delta bias and delta uncertainty computations


%% ==> choice predictivity and AIC analysis