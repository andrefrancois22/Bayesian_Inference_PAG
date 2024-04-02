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

% ==> delta bias and uncertainty metrics
dbias_lo_ctr = nan(29,1);
dpu_lo_ctr   = nan(29,1);

% ==> PF curve fits
PFs_predPF_dynr_lo = {};
PFs_predPF_dynr_hi = {};

figure(4); set(gcf,'color','white'); set(gcf,'Position',[91 109 1160 853]);

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
    
    % ==> choose stimulus
    stim_idx = true(length(ori),1);
%     stim_idx = (ori == 0);;
%     stim_idx = ( abs(ori) >= abs(or(2)) );
    
    % => indexing variables - dynamic range
    lo_dynr = (dynr(:,1) <  median(dynr(stim_idx,1)));    
    hi_dynr = (dynr(:,1) >= median(dynr(stim_idx,1)));

    % ==> cw context (congruent choices)
    % ==> cw context & lo contrast, and lo dynamic range
    cw_lcr_ldr_cng = (ctx == 1) & (ctr == min(ctr)) & (lo_dynr) & (cho == 1) & stim_idx;    
    % ==> cw context & lo contrast, and hi dynamic range
    cw_lcr_hdr_cng = (ctx == 1) & (ctr == min(ctr)) & (hi_dynr) & (cho == 1) & stim_idx;    
    % ==> ccw context (congruent choices)
    % ==> ccw context & lo contrast, and lo dynamic range
    ccw_lcr_ldr_cng = (ctx == -1) & (ctr == min(ctr)) & (lo_dynr) & (cho == -1) & stim_idx;    
    % ==> ccw context & lo contrast, and hi dynamic range
    ccw_lcr_hdr_cng = (ctx == -1) & (ctr == min(ctr)) & (hi_dynr) & (cho == -1) & stim_idx;      
        
    
    % ==> cw context (incongruent choices)
    % ==> cw context & lo contrast, and lo dynamic range
    cw_lcr_ldr_icg = (ctx == 1) & (ctr == min(ctr)) & (lo_dynr) & (cho == -1) & stim_idx;    
    % ==> cw context & lo contrast, and hi dynamic range
    cw_lcr_hdr_icg = (ctx == 1) & (ctr == min(ctr)) & (hi_dynr) & (cho == -1) & stim_idx;            
    % ==> ccw context (incongruent choices)
    % ==> ccw context & lo contrast, and lo dynamic range
    ccw_lcr_ldr_icg = (ctx == -1) & (ctr == min(ctr)) & (lo_dynr) & (cho == 1) & stim_idx;    
    % ==> ccw context & lo contrast, and hi dynamic range
    ccw_lcr_hdr_icg = (ctx == -1) & (ctr == min(ctr)) & (hi_dynr) & (cho == 1) & stim_idx;       
    
    % ==> counts
    % ==> lo contrast stimuli
    cho_cw_ct_cw_cr_lo_dynr_lo   = nan(1,7);
    cho_cw_ct_ccw_cr_lo_dynr_lo  = nan(1,7);
    
    cho_ccw_ct_cw_cr_lo_dynr_lo  = nan(1,7);
    cho_ccw_ct_ccw_cr_lo_dynr_lo = nan(1,7);
    
    cho_cw_ct_cw_cr_lo_dynr_hi   = nan(1,7);
    cho_cw_ct_ccw_cr_lo_dynr_hi  = nan(1,7);
    
    cho_ccw_ct_cw_cr_lo_dynr_hi  = nan(1,7);
    cho_ccw_ct_ccw_cr_lo_dynr_hi = nan(1,7);
    
    % ==> counts
    % ==> high contrast stimuli
    cho_cw_ct_cw_cr_hi_dynr_lo   = nan(1,7);
    cho_cw_ct_ccw_cr_hi_dynr_lo  = nan(1,7);
    
    cho_ccw_ct_cw_cr_hi_dynr_lo  = nan(1,7);
    cho_ccw_ct_ccw_cr_hi_dynr_lo = nan(1,7);
    
    cho_cw_ct_cw_cr_hi_dynr_hi   = nan(1,7);
    cho_cw_ct_ccw_cr_hi_dynr_hi  = nan(1,7);
    
    cho_ccw_ct_cw_cr_hi_dynr_hi  = nan(1,7);
    cho_ccw_ct_ccw_cr_hi_dynr_hi = nan(1,7);    
    
    % ==> compute counts for each stimulus orientation
    for theta = 1:7
        
        % ==> low contrast trials only
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

    
        % ==> high contrast trials only
    
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

    % ==> low contrast
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
    
    % ==> high contrast
    
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
    dbias_lo_ctr(iSfln) = (paramEst_M1_dynr_lo(5) - paramEst_M1_dynr_lo(4)) - (paramEst_M1_dynr_hi(5) - paramEst_M1_dynr_hi(4));
    %paramEst_M1_dynr_lo(5) - paramEst_M1_dynr_hi(5);    
    % ==> delta perceptual uncertainty
    dpu_lo_ctr(iSfln)   = paramEst_M1_dynr_lo(3) - paramEst_M1_dynr_hi(3);
    
    % ==> save PF curve fits
    PFs_predPF_dynr_lo{iSfln} = predPF_dynr_lo;
    PFs_predPF_dynr_hi{iSfln} = predPF_dynr_hi;    
    
end


% ==> DV choice predictivity
dvCatPerf = load([drc,'dvCatPerf.mat']);
dvCatPerf = dvCatPerf.dvCatPerf';

% ==> AIC Delta
aic_d = load([drc,'aic_d.mat']);
aic_d = aic_d.aic_d';

% ==> log-likelihood ratios
llr = load([drc,'llr.mat']);
llr = llr.llr';

% ==> choose populations to plot
ppl = 2:29;
% ppl(ppl==19) = [];

% ==> choice predictivity and delta bias + delta uncertainty
[rb,pb] = corr(dvCatPerf(ppl), dbias_lo_ctr(ppl));
[rp,pp] = corr(dvCatPerf(ppl), dpu_lo_ctr(ppl));

% ==> Use AIC computation from F2 instead of choice predictivity
[rb_aic,pb_aic] = corr(aic_d(ppl), dbias_lo_ctr(ppl));
[rp_aic,pp_aic] = corr(aic_d(ppl), dpu_lo_ctr(ppl));

% ==> with LLR from F2
[rb_llr,pb_llr] = corr(llr(ppl), dbias_lo_ctr(ppl));
[rp_llr,pp_llr] = corr(llr(ppl), dpu_lo_ctr(ppl));


% ==> pearson correlation between delta bias and delta uncertainty
% ==> ommitting session 1 here!
[r_bp,p_bp] = corr(dbias_lo_ctr(ppl),dpu_lo_ctr(ppl));

figure(5); set(gcf,'color','white'); set(gcf,'Position',[723 4 1104 958]);
subplot(3,3,1);
scatter(dvCatPerf(ppl), dbias_lo_ctr(ppl),30,'ko','filled'); hold on; hold all;
text(dvCatPerf(ppl)+0.005, dbias_lo_ctr(ppl)-0.005, string([ppl]),'FontSize',6); 
title(['r = ',num2str(rb), ', p = ',num2str(pb)]);
xlabel('Choice Predictivity');
ylabel('\Delta bias');

subplot(3,3,2);
scatter(dvCatPerf(ppl), dpu_lo_ctr(ppl),30,'ko','filled'); hold on; hold all;
text(dvCatPerf(ppl)+0.005, dpu_lo_ctr(ppl)-0.005,string([ppl]),'FontSize',6);
title(['r = ',num2str(rp), ', p = ',num2str(pp)]);
xlabel('Choice Predictivity');
ylabel('\Delta perceptual uncertainty (slope)');

subplot(3,3,3);
scatter(dbias_lo_ctr(ppl), dpu_lo_ctr(ppl),30,'ko','filled'); hold on; hold all;
text(dbias_lo_ctr(ppl)+0.005,dpu_lo_ctr(ppl)-0.005,string([ppl]),'FontSize',6);
title(['r = ',num2str(r_bp), ', p = ',num2str(p_bp)]);
xlabel('\Delta bias')
ylabel('\Delta perceptual uncertainty (slope)')

subplot(3,3,4);
scatter(aic_d(ppl), dbias_lo_ctr(ppl),30,'ko','filled'); hold on; hold all;
text(aic_d(ppl)+0.005, dbias_lo_ctr(ppl)-0.005, string([ppl]),'FontSize',6); 
title(['r = ',num2str(rb_aic), ', p = ',num2str(pb_aic)]);
xlabel('AIC \Delta');
ylabel('\Delta bias');

subplot(3,3,5);
scatter(aic_d(ppl), dpu_lo_ctr(ppl),30,'ko','filled'); hold on; hold all;
text(aic_d(ppl)+0.005, dpu_lo_ctr(ppl)-0.005,string([ppl]),'FontSize',6);
title(['r = ',num2str(rp_aic), ', p = ',num2str(pp_aic)]);
xlabel('AIC \Delta');
ylabel('\Delta perceptual uncertainty (slope)');

subplot(3,3,7);
scatter(llr(ppl), dbias_lo_ctr(ppl),30,'ko','filled'); hold on; hold all;
text(llr(ppl)+0.005, dbias_lo_ctr(ppl)-0.005, string([ppl]),'FontSize',6); 
title(['r = ',num2str(rb_llr), ', p = ',num2str(pb_llr)]);
xlabel('Log-likelihood ratio');
ylabel('\Delta bias');

subplot(3,3,8);
scatter(llr(ppl), dpu_lo_ctr(ppl),30,'ko','filled'); hold on; hold all;
text(llr(ppl)+0.005, dpu_lo_ctr(ppl)-0.005,string([ppl]),'FontSize',6);
title(['r = ',num2str(rp_llr), ', p = ',num2str(pp_llr)]);
xlabel('Log-likelihood ratio');
ylabel('\Delta perceptual uncertainty (slope)');

%% ==>

PFs_predPF_dynr_lo_mu = nan(29,2,7);

cd = colormap('parula'); % take your pick (doc colormap)
cd = interp1(linspace(min(aic_d(ppl)),max(aic_d),length(cd)),cd,aic_d(ppl)); % map color to y values
cd = uint8(cd'*255); % need a 4xN uint8 array
cd(4,:) = 255; % last column is transparency


figure(6); set(gcf,'color','white');
% subplot(1,2,1);
for iS = ppl
    plot([-3,-2,-1,0,1,2,3], PFs_predPF_dynr_lo{iS}(1,:) - PFs_predPF_dynr_lo{iS}(2,:), '-', 'Color',cd(:,iS-1), 'LineWidth', 1); hold on; hold all;
    xlabel('Orientation')
    ylabel('\Delta bias')
    PFs_predPF_dynr_lo_mu(iS,1,:) = PFs_predPF_dynr_lo{iS}(1,:);
    PFs_predPF_dynr_lo_mu(iS,2,:) = PFs_predPF_dynr_lo{iS}(2,:);
end
colorbar();

% ==> average
mu = squeeze(nanmean(PFs_predPF_dynr_lo_mu,1));
y = mu(1,:)-mu(2,:);
h = plot(or,y,'k.-','LineWidth',2);

dfs = squeeze((PFs_predPF_dynr_lo_mu(ppl,1,:) - PFs_predPF_dynr_lo_mu(ppl,2,:)));



% figure();
% plot(or, ncw_dynr_lo./sum(ncw_dynr_lo,1));
% hold on; hold all;
% plot(or, nccw_dynr_lo./sum(nccw_dynr_lo,1));
% 
% ncw_dynr_hi./sum(ncw_dynr_hi,1)
% nccw_dynr_hi./sum(nccw_dynr_hi,1)



