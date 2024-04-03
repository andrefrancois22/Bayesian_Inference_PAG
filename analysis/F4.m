clear all; close all; clc;
% ==> directories
dataPath     = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data';
drc = '../data/';
% ==> directory containing PF curve fit functions
pfc_functions_dr = '/home/thomas/Desktop/Bayesian_Inference_PAG/simulation/pfc_functions/';
addpath(pfc_functions_dr)

%% ==> compute DV dynamic range medians for each stimulus type (29, 7 x 2 x 2) e.g. population x orientation x context by contrast

% ==> DV dynamic range medians for each stimulus
DV_dynr_meds = nan(29,7,2,2);

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

    % ==> load session data
    if iS > 9
        load([dataPath,'/dataSet_',num2str(iS),'.mat']);
    elseif iS <= 9
        load([dataPath,'/dataSet_0',num2str(iS),'.mat']);
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
    % dynr = max([max(dvs, [], 2) - dvs(:,1), -(min(dvs, [], 2) - dvs(:,1))], [], 2);

    for oi = 1:length(or)
        for xi = 1:length(cx)
            for ci = 1:length(cr)                
                % ==> contrsuct stimulus index
                I = (ori == or(oi)) & (ctx == cx(xi)) & (ctr == cr(ci));                
                % ==> get DV trajectories indexed by stimulus
                dv = dvs(I,:);
                % ==> compute dynamic range for each trial for DV traj. of
                % that stimulus
                dynr = max([max(dv, [], 2) - dv(:,1), -(min(dv, [], 2) - dv(:,1))], [], 2);
                % ==> compute the median
                md = median(dynr);
                % store the dyn. range median for that stimulus
                DV_dynr_meds(iS, oi, xi, ci) = md;
            end
        end
    end
    % ==> update
    fprintf('computed all DV dynamic range medians for population %d of %d...\n',iS,29)
end
fprintf('Done...\n')

%% ==> dynamic range analysis

% ==> dynamic range function
dynf = @(x) max([max(x, [], 2) - x(:,1), -(min(x, [], 2) - x(:,1))], [], 2);

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
db_lcr = nan(29,1);
dp_lcr   = nan(29,1);
% ==> high contrast cases
db_hcr = nan(29,1);
dp_hcr   = nan(29,1);

% ==> PF curve fits
PFs_predPF_dynr_lo = cell(29,1);
PFs_predPF_dynr_hi = cell(29,1);

% ==> PF curve fits
PFs_predPF_dynr_lo_hcr = cell(29,1);
PFs_predPF_dynr_hi_hcr = cell(29,1);

figure(3); set(gcf,'color','white'); set(gcf,'Position',[91 109 1160 853]);
figure(4); set(gcf,'color','white'); set(gcf,'Position',[91 109 1160 853]);

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

    % ==> load session data
    if iS > 9
        load([dataPath,'/dataSet_',num2str(iS),'.mat']);
    elseif iS <= 9
        load([dataPath,'/dataSet_0',num2str(iS),'.mat']);
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
    dynr = dynf(dvs);        
    % => indexing variables - dynamic range
    ldyn = (dynr(:,1) <  median(dynr));    
    hdyn = (dynr(:,1) >= median(dynr));

    % ==> cw context (congruent choices)
    % ==> cw context & lo contrast, and lo dynamic range
    cw_lcr_ldr_cng = (ctx == 1) & (ctr == min(ctr)) & (ldyn) & (cho == 1); % ==> DV_dynr_meds(iS,th,2,1)
    % ==> cw context & lo contrast, and hi dynamic range
    cw_lcr_hdr_cng = (ctx == 1) & (ctr == min(ctr)) & (hdyn) & (cho == 1); % ==> DV_dynr_meds(iS,th,2,1)    
    % ==> ccw context (congruent choices)
    % ==> ccw context & lo contrast, and lo dynamic range
    ccw_lcr_ldr_cng = (ctx == -1) & (ctr == min(ctr)) & (ldyn) & (cho == -1); %  
    % ==> ccw context & lo contrast, and hi dynamic range
    ccw_lcr_hdr_cng = (ctx == -1) & (ctr == min(ctr)) & (hdyn) & (cho == -1); %      
    % ==> cw context (incongruent choices)
    % ==> cw context & lo contrast, and lo dynamic range
    cw_lcr_ldr_icg = (ctx == 1) & (ctr == min(ctr)) & (ldyn) & (cho == -1); %   
    % ==> cw context & lo contrast, and hi dynamic range
    cw_lcr_hdr_icg = (ctx == 1) & (ctr == min(ctr)) & (hdyn) & (cho == -1); %         
    % ==> ccw context (incongruent choices)
    % ==> ccw context & lo contrast, and lo dynamic range
    ccw_lcr_ldr_icg = (ctx == -1) & (ctr == min(ctr)) & (ldyn) & (cho == 1); %  
    % ==> ccw context & lo contrast, and hi dynamic range
    ccw_lcr_hdr_icg = (ctx == -1) & (ctr == min(ctr)) & (hdyn) & (cho == 1); %      
    
    
    % ==> cw context & lo contrast, and lo dynamic range
    cw_hcr_ldr_cng = (ctx == 1) & (ctr == max(ctr)) & (ldyn) & (cho == 1); %    
    % ==> cw context & lo contrast, and hi dynamic range
    cw_hcr_hdr_cng = (ctx == 1) & (ctr == max(ctr)) & (hdyn) & (cho == 1); %   
    % ==> ccw context (congruent choices)
    % ==> ccw context & lo contrast, and lo dynamic range
    ccw_hcr_ldr_cng = (ctx == -1) & (ctr == max(ctr)) & (ldyn) & (cho == -1); %    
    % ==> ccw context & lo contrast, and hi dynamic range
    ccw_hcr_hdr_cng = (ctx == -1) & (ctr == max(ctr)) & (hdyn) & (cho == -1); %     
    % ==> cw context (incongruent choices)
    % ==> cw context & lo contrast, and lo dynamic range
    cw_hcr_ldr_icg = (ctx == 1) & (ctr == max(ctr)) & (ldyn) & (cho == -1); %    
    % ==> cw context & lo contrast, and hi dynamic range
    cw_hcr_hdr_icg = (ctx == 1) & (ctr == max(ctr)) & (hdyn) & (cho == -1); %           
    % ==> ccw context (incongruent choices)
    % ==> ccw context & lo contrast, and lo dynamic range
    ccw_hcr_ldr_icg = (ctx == -1) & (ctr == max(ctr)) & (ldyn) & (cho == 1); %   
    % ==> ccw context & lo contrast, and hi dynamic range
    ccw_hcr_hdr_icg = (ctx == -1) & (ctr == max(ctr)) & (hdyn) & (cho == 1); %       
                       
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
    for th = 1:length(or)
        
        % ==> low contrast trials only ==================================================================================================================
        % ==> low dynamic range split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % ==> cw congruent choice, lo ctr, lo dynr
        
%         % ==> dynamic range for session
%         dynr = dynf(dvs);    
%         % => indexing variables - dynamic range
%         ldyn = (dynr(:,1) <  median(dynr)); hdyn = (dynr(:,1) >= median(dynr));
        
        cho_cw_ct_cw_cr_lo_dynr_lo(th)   = sum(cw_lcr_ldr_cng   & ori==or(th));
        % ==> cw incongruent choice, lo ctr, lo dynr
        cho_cw_ct_ccw_cr_lo_dynr_lo(th)  = sum(ccw_lcr_ldr_icg  & ori==or(th));      
        % ==> ccw incongruent choice, lo ctr, lo dynr
        cho_ccw_ct_cw_cr_lo_dynr_lo(th)  = sum(cw_lcr_ldr_icg   & ori==or(th));
        % ==> ccw congruent choice, lo ctr, lo dynr
        cho_ccw_ct_ccw_cr_lo_dynr_lo(th) = sum(ccw_lcr_ldr_cng  & ori==or(th));      
        
        % ==> high dynamic range split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         % ==> cw congruent choice, lo ctr, lo dynr
        cho_cw_ct_cw_cr_lo_dynr_hi(th)   = sum(cw_lcr_hdr_cng   & ori==or(th));
        % ==> cw incongruent choice, lo ctr, lo dynr
        cho_cw_ct_ccw_cr_lo_dynr_hi(th)  = sum(ccw_lcr_hdr_icg  & ori==or(th));      
        % ==> ccw incongruent choice, lo ctr, lo dynr
        cho_ccw_ct_cw_cr_lo_dynr_hi(th)  = sum(cw_lcr_hdr_icg   & ori==or(th));
        % ==> ccw congruent choice, lo ctr, lo dynr
        cho_ccw_ct_ccw_cr_lo_dynr_hi(th) = sum(ccw_lcr_hdr_cng  & ori==or(th));
        % ================================================================================================================================================

    
        % ==> high contrast trials only =================================================================================================================
        % ==> low dynamic range split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % ==> cw congruent choice, lo ctr, lo dynr
        cho_cw_ct_cw_cr_hi_dynr_lo(th)   = sum(cw_hcr_ldr_cng   & ori==or(th));
        % ==> cw incongruent choice, lo ctr, lo dynr
        cho_cw_ct_ccw_cr_hi_dynr_lo(th)  = sum(ccw_hcr_ldr_icg  & ori==or(th));      
        % ==> ccw incongruent choice, lo ctr, lo dynr
        cho_ccw_ct_cw_cr_hi_dynr_lo(th)  = sum(cw_hcr_ldr_icg   & ori==or(th));
        % ==> ccw congruent choice, lo ctr, lo dynr
        cho_ccw_ct_ccw_cr_hi_dynr_lo(th) = sum(ccw_hcr_ldr_cng  & ori==or(th));
        
        % ==> high dynamic range split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         % ==> cw congruent choice, lo ctr, lo dynr
        cho_cw_ct_cw_cr_hi_dynr_hi(th)   = sum(cw_hcr_hdr_cng   & ori==or(th));
        % ==> cw incongruent choice, lo ctr, lo dynr
        cho_cw_ct_ccw_cr_hi_dynr_hi(th)  = sum(ccw_hcr_hdr_icg  & ori==or(th));        
        % ==> ccw incongruent choice, lo ctr, lo dynr
        cho_ccw_ct_cw_cr_hi_dynr_hi(th)  = sum(cw_hcr_hdr_icg   & ori==or(th));
        % ==> ccw congruent choice, lo ctr, lo dynr
        cho_ccw_ct_ccw_cr_hi_dynr_hi(th) = sum(ccw_hcr_hdr_cng  & ori==or(th));    
        % ================================================================================================================================================
    
    end
    
    % ==> PF curve fit
    fprintf('computing psychometric functions...\n')    
    % ==> low dynamic range cases (low contrast stimulus trials)
    [PF_ldyn_lcr, p_ldyn_lcr] = PF_fit_fun(cho_ccw_ct_cw_cr_lo_dynr_lo, cho_ccw_ct_ccw_cr_lo_dynr_lo, cho_cw_ct_cw_cr_lo_dynr_lo, cho_cw_ct_ccw_cr_lo_dynr_lo, or, startVec_M1, LB_M1, UB_M1, options);
    % ==> high dynamic range cases (low contrast stimulus trials)
    [PF_hdyn_lcr, p_hdyn_lcr] = PF_fit_fun(cho_ccw_ct_cw_cr_lo_dynr_hi, cho_ccw_ct_ccw_cr_lo_dynr_hi, cho_cw_ct_cw_cr_lo_dynr_hi, cho_cw_ct_ccw_cr_lo_dynr_hi, or, startVec_M1, LB_M1, UB_M1, options);
    
    % ==> low dynamic range cases (high contrast stimulus trials)
    [PF_ldyn_hcr, p_ldyn_hcr] = PF_fit_fun(cho_ccw_ct_cw_cr_hi_dynr_lo, cho_ccw_ct_ccw_cr_hi_dynr_lo, cho_cw_ct_cw_cr_hi_dynr_lo, cho_cw_ct_ccw_cr_hi_dynr_lo, or, startVec_M1, LB_M1, UB_M1, options);
    % ==> high dynamic range cases (high contrast stimulus trials)        
    [PF_hdyn_hcr, p_hdyn_hcr] = PF_fit_fun(cho_ccw_ct_cw_cr_hi_dynr_hi, cho_ccw_ct_ccw_cr_hi_dynr_hi, cho_cw_ct_cw_cr_hi_dynr_hi, cho_cw_ct_ccw_cr_hi_dynr_hi, or, startVec_M1, LB_M1, UB_M1, options);
    
    % ==> plots
    figure(3);
    subplot(5,6,iS); title(['Session: ', num2str(iS)]);
    plot(or,PF_ldyn_lcr(1,:),'b.-')
    hold on; hold all;
    plot(or,PF_ldyn_lcr(2,:),'r.-')
    hold on; hold all;
    plot(or,PF_hdyn_lcr(1,:),'b+:')
    hold on; hold all;
    plot(or,PF_hdyn_lcr(2,:),'r+:')    
    xlabel('stimulus orientation')
    ylabel('Proportion cw PF fit')
    drawnow;    
    % ==> plots
    figure(4);
    subplot(5,6,iS); title(['Session: ', num2str(iS)]);
    plot(or,PF_ldyn_hcr(1,:),'g.-')
    hold on; hold all;
    plot(or,PF_ldyn_hcr(2,:),'m.-')
    hold on; hold all;
    plot(or,PF_hdyn_hcr(1,:),'g+:')
    hold on; hold all;
    plot(or,PF_hdyn_hcr(2,:),'m+:')    
    xlabel('stimulus orientation')
    ylabel('Proportion cw PF fit')
    drawnow;    
    
    % ==> delta bias (this is a different computation from original)
    db_lcr(iS) = (p_ldyn_lcr(5) - p_ldyn_lcr(4)) - (p_hdyn_lcr(5) - p_hdyn_lcr(4));
    %paramEst_M1_dynr_lo(5) - paramEst_M1_dynr_hi(5);    
    % ==> delta perceptual uncertainty
    dp_lcr(iS)   = p_ldyn_lcr(3) - p_hdyn_lcr(3);
    
    % ==> delta bias (this is a different computation from original)
    db_hcr(iS) = (p_ldyn_hcr(5) - p_ldyn_hcr(4)) - (p_hdyn_hcr(5) - p_hdyn_hcr(4));
    %paramEst_M1_dynr_lo_hcr(5) - paramEst_M1_dynr_hi)_hcr(5);    
    % ==> delta perceptual uncertainty
    dp_hcr(iS)   = p_ldyn_hcr(3) - p_hdyn_hcr(3);    
    
    % ==> save PF curve fits
    PFs_predPF_dynr_lo{iS} = PF_ldyn_lcr;
    PFs_predPF_dynr_hi{iS} = PF_hdyn_lcr;  
    % ==> save PF curve fits
    PFs_predPF_dynr_lo_hcr{iS} = PF_ldyn_hcr;
    PFs_predPF_dynr_hi_hcr{iS} = PF_hdyn_hcr;      
    
end

%% ==> DV choice predictivity

dvCatPerf = load([drc,'dvCatPerf.mat']);
dvCatPerf = dvCatPerf.dvCatPerf';

% dbias_lo_ctr(dbias_lo_ctr > 4) = 4;
% dbias_lo_ctr_hcr(dbias_lo_ctr_hcr > 4) = 4;

% ==> AIC Delta
aic_d = load([drc,'aic_d.mat']);
aic_d = aic_d.aic_d';

% ==> log-likelihood ratios
llr = load([drc,'llr.mat']);
llr = llr.llr';

% ==> choose populations to plot
% ppl = 1:13 %14:29 
ppl = 2:29;
% ppl(ppl==19) = [];

% ==> choice predictivity and delta bias + delta uncertainty
[rb,pb] = corr(dvCatPerf(ppl), db_lcr(ppl));
[rp,pp] = corr(dvCatPerf(ppl), dp_lcr(ppl));

% ==> Use AIC computation from F2 instead of choice predictivity
[rb_aic,pb_aic] = corr(aic_d(ppl), db_lcr(ppl));
[rp_aic,pp_aic] = corr(aic_d(ppl), dp_lcr(ppl));

% ==> with LLR from F2
[rb_llr,pb_llr] = corr(llr(ppl), db_lcr(ppl));
[rp_llr,pp_llr] = corr(llr(ppl), dp_lcr(ppl));


% ==> pearson correlation between delta bias and delta uncertainty
% ==> ommitting session 1 here!
[r_bp,p_bp] = corr(db_lcr(ppl),dp_lcr(ppl));


figure(5); set(gcf,'color','white'); set(gcf,'Position',[723 4 1104 958]);
subplot(3,3,1);
scatter(dvCatPerf(ppl), db_lcr(ppl),30,'ko','filled'); hold on; hold all;
text(dvCatPerf(ppl)+0.005, db_lcr(ppl)-0.005, string([ppl]),'FontSize',6); 
title(['r = ',num2str(rb), ', p = ',num2str(pb)]);
xlabel('Choice Predictivity');
ylabel('\Delta bias');
subplot(3,3,2);
scatter(dvCatPerf(ppl), dp_lcr(ppl),30,'ko','filled'); hold on; hold all;
text(dvCatPerf(ppl)+0.005, dp_lcr(ppl)-0.005,string([ppl]),'FontSize',6);
title(['r = ',num2str(rp), ', p = ',num2str(pp)]);
xlabel('Choice Predictivity');
ylabel('\Delta perceptual uncertainty (slope)');
subplot(3,3,3);
scatter(db_lcr(ppl), dp_lcr(ppl),30,'ko','filled'); hold on; hold all;
text(db_lcr(ppl)+0.005,dp_lcr(ppl)-0.005,string([ppl]),'FontSize',6);
title(['r = ',num2str(r_bp), ', p = ',num2str(p_bp)]);
xlabel('\Delta bias')
ylabel('\Delta perceptual uncertainty (slope)')
subplot(3,3,4);
scatter(aic_d(ppl), db_lcr(ppl),30,'ko','filled'); hold on; hold all;
text(aic_d(ppl)+0.005, db_lcr(ppl)-0.005, string([ppl]),'FontSize',6); 
title(['r = ',num2str(rb_aic), ', p = ',num2str(pb_aic)]);
xlabel('AIC \Delta');
ylabel('\Delta bias');
subplot(3,3,5);
scatter(aic_d(ppl), dp_lcr(ppl),30,'ko','filled'); hold on; hold all;
text(aic_d(ppl)+0.005, dp_lcr(ppl)-0.005,string([ppl]),'FontSize',6);
title(['r = ',num2str(rp_aic), ', p = ',num2str(pp_aic)]);
xlabel('AIC \Delta');
ylabel('\Delta perceptual uncertainty (slope)');
subplot(3,3,7);
scatter(llr(ppl), db_lcr(ppl),30,'ko','filled'); hold on; hold all;
text(llr(ppl)+0.005, db_lcr(ppl)-0.005, string([ppl]),'FontSize',6); 
title(['r = ',num2str(rb_llr), ', p = ',num2str(pb_llr)]);
xlabel('Log-likelihood ratio');
ylabel('\Delta bias');
subplot(3,3,8);
scatter(llr(ppl), dp_lcr(ppl),30,'ko','filled'); hold on; hold all;
text(llr(ppl)+0.005, dp_lcr(ppl)-0.005,string([ppl]),'FontSize',6);
title(['r = ',num2str(rp_llr), ', p = ',num2str(pp_llr)]);
xlabel('Log-likelihood ratio');
ylabel('\Delta perceptual uncertainty (slope)');

%% ==> high contrast trials

% ==> choose populations to plot
% ppl = 1:13% 14:29 
ppl = 2:29;

% ==> choice predictivity and delta bias + delta uncertainty
[rb,pb] = corr(dvCatPerf(ppl), db_hcr(ppl));
[rp,pp] = corr(dvCatPerf(ppl), dp_hcr(ppl));

% ==> Use AIC computation from F2 instead of choice predictivity
[rb_aic,pb_aic] = corr(aic_d(ppl), db_hcr(ppl));
[rp_aic,pp_aic] = corr(aic_d(ppl), dp_hcr(ppl));

% ==> with LLR from F2
[rb_llr,pb_llr] = corr(llr(ppl), db_hcr(ppl));
[rp_llr,pp_llr] = corr(llr(ppl), dp_hcr(ppl));


% ==> pearson correlation between delta bias and delta uncertainty
% ==> ommitting session 1 here!
[r_bp,p_bp] = corr(db_hcr(ppl),dp_hcr(ppl));

figure(6); set(gcf,'color','white'); set(gcf,'Position',[723 4 1104 958]);
subplot(3,3,1);
scatter(dvCatPerf(ppl), db_hcr(ppl),30,'mo','filled'); hold on; hold all;
text(dvCatPerf(ppl)+0.005, db_hcr(ppl)-0.005, string([ppl]),'FontSize',6); 
title(['r = ',num2str(rb), ', p = ',num2str(pb)]);
xlabel('Choice Predictivity');
ylabel('\Delta bias');
subplot(3,3,2);
scatter(dvCatPerf(ppl), dp_hcr(ppl),30,'mo','filled'); hold on; hold all;
text(dvCatPerf(ppl)+0.005, dp_hcr(ppl)-0.005,string([ppl]),'FontSize',6);
title(['r = ',num2str(rp), ', p = ',num2str(pp)]);
xlabel('Choice Predictivity');
ylabel('\Delta perceptual uncertainty (slope)');
subplot(3,3,3);
scatter(db_hcr(ppl), dp_hcr(ppl),30,'mo','filled'); hold on; hold all;
text(db_hcr(ppl)+0.005,dp_hcr(ppl)-0.005,string([ppl]),'FontSize',6);
title(['r = ',num2str(r_bp), ', p = ',num2str(p_bp)]);
xlabel('\Delta bias')
ylabel('\Delta perceptual uncertainty (slope)')
subplot(3,3,4);
scatter(aic_d(ppl), db_hcr(ppl),30,'mo','filled'); hold on; hold all;
text(aic_d(ppl)+0.005, db_hcr(ppl)-0.005, string([ppl]),'FontSize',6); 
title(['r = ',num2str(rb_aic), ', p = ',num2str(pb_aic)]);
xlabel('AIC \Delta');
ylabel('\Delta bias');
subplot(3,3,5);
scatter(aic_d(ppl), dp_hcr(ppl),30,'mo','filled'); hold on; hold all;
text(aic_d(ppl)+0.005, dp_hcr(ppl)-0.005,string([ppl]),'FontSize',6);
title(['r = ',num2str(rp_aic), ', p = ',num2str(pp_aic)]);
xlabel('AIC \Delta');
ylabel('\Delta perceptual uncertainty (slope)');
subplot(3,3,7);
scatter(llr(ppl), db_hcr(ppl),30,'mo','filled'); hold on; hold all;
text(llr(ppl)+0.005, db_hcr(ppl)-0.005, string([ppl]),'FontSize',6); 
title(['r = ',num2str(rb_llr), ', p = ',num2str(pb_llr)]);
xlabel('Log-likelihood ratio');
ylabel('\Delta bias');
subplot(3,3,8);
scatter(llr(ppl), dp_hcr(ppl),30,'mo','filled'); hold on; hold all;
text(llr(ppl)+0.005, dp_hcr(ppl)-0.005,string([ppl]),'FontSize',6);
title(['r = ',num2str(rp_llr), ', p = ',num2str(pp_llr)]);
xlabel('Log-likelihood ratio');
ylabel('\Delta perceptual uncertainty (slope)');

%% ==>

% PFs_predPF_dynr_lo_mu = nan(29,2,7);
% 
% cd = colormap('parula'); % take your pick (doc colormap)
% cd = interp1(linspace(min(aic_d(ppl)),max(aic_d),length(cd)),cd,aic_d(ppl)); % map color to y values
% cd = uint8(cd'*255); % need a 4xN uint8 array
% cd(4,:) = 255; % last column is transparency
% 
% 
% figure(6); set(gcf,'color','white');
% % subplot(1,2,1);
% for iS = ppl
%     plot([-3,-2,-1,0,1,2,3], PFs_predPF_dynr_lo{iS}(1,:) - PFs_predPF_dynr_lo{iS}(2,:), '-', 'Color',cd(:,iS-1), 'LineWidth', 1); hold on; hold all;
%     xlabel('Orientation')
%     ylabel('\Delta bias')
%     PFs_predPF_dynr_lo_mu(iS,1,:) = PFs_predPF_dynr_lo{iS}(1,:);
%     PFs_predPF_dynr_lo_mu(iS,2,:) = PFs_predPF_dynr_lo{iS}(2,:);
% end
% colorbar();
% 
% % ==> average
% mu = squeeze(nanmean(PFs_predPF_dynr_lo_mu,1));
% y = mu(1,:)-mu(2,:);
% h = plot(or,y,'k.-','LineWidth',2);
% 
% dfs = squeeze((PFs_predPF_dynr_lo_mu(ppl,1,:) - PFs_predPF_dynr_lo_mu(ppl,2,:)));

% figure();
% plot(or, ncw_dynr_lo./sum(ncw_dynr_lo,1));
% hold on; hold all;
% plot(or, nccw_dynr_lo./sum(nccw_dynr_lo,1));
% 
% ncw_dynr_hi./sum(ncw_dynr_hi,1)
% nccw_dynr_hi./sum(nccw_dynr_hi,1)



