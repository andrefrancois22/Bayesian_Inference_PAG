clear all; close all; clc;
% ==> directories
dataPath     = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data';
drc = '../data/';
% ==> directory containing PF curve fit functions
pfc_functions_dr = '/home/thomas/Desktop/Bayesian_Inference_PAG/simulation/F1/pfc_functions/';
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
close all; clc;

% ==> estimate dynamic range per population?
dyn_type = 'popu';
% ==> estimate dynamic range per stimulus?
% dyn_type = 'stim';

DV_dynr_meds_check = nan(29,7,2,2);

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
LB_M1(3,1)  = 0.1;                       UB_M1(3,1) = 10;        % perceptual uncertainty
LB_M1(4,1)  = -5;                        UB_M1(4,1) = 5;        % decision criterion (-10->10 range originally)
LB_M1(5,1)  = -5;                        UB_M1(5,1) = 5;        % decision criterion
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

figure(3); set(gcf,'color','white'); set(gcf,'Position',[86 282 921 668]);
% figure(4); set(gcf,'color','white'); set(gcf,'Position',[86 282 921 668]);
% figure(5); set(gcf,'color','white'); set(gcf,'Position',[86 282 921 668]);
% figure(6); set(gcf,'color','white'); set(gcf,'Position',[86 282 921 668]);


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
    
    if strcmp(dyn_type,'popu')
        % ==> dynamic range for session
        dynr = dynf(dvs); ldyn = (dynr(:,1) <  median(dynr)); hdyn = (dynr(:,1) >= median(dynr));
    end
    
    % ==> cw context (congruent choices)
    % ==> cw context & lo contrast
    cw_lcr_cng = (ctx == 1) & (ctr == min(ctr))   & (cho == 1); 
    % ==> ccw context & lo contrast
    ccw_lcr_cng = (ctx == -1) & (ctr == min(ctr)) & (cho == -1); %     
    % ==> cw context (incongruent choices)
    % ==> cw context & lo contrast
    cw_lcr_icg = (ctx == 1) & (ctr == min(ctr))   & (cho == -1); %          
    % ==> ccw context (incongruent choices)
    % ==> ccw context & lo contrast
    ccw_lcr_icg = (ctx == -1) & (ctr == min(ctr)) & (cho == 1); %   
        
    % ==> cw context & hi contrast
    cw_hcr_cng = (ctx == 1) & (ctr == max(ctr))   & (cho == 1); %    
    % ==> ccw context (congruent choices)
    % ==> ccw context & hi contrast
    ccw_hcr_cng = (ctx == -1) & (ctr == max(ctr)) & (cho == -1); %     
    % ==> cw context (incongruent choices)
    % ==> cw context & hi contrast
    cw_hcr_icg = (ctx == 1) & (ctr == max(ctr))   & (cho == -1); %             
    % ==> ccw context (incongruent choices)
    % ==> ccw context & hi contrast
    ccw_hcr_icg = (ctx == -1) & (ctr == max(ctr)) & (cho == 1); %      
                       
    % ==> counts
    % ==> lo contrast stimuli (==> <choice>_<context>_<contrast>_<dynamic range>)
    cw_cw_lo_lo   = nan(1,7); 
    cw_ccw_lo_lo  = nan(1,7);    
    ccw_cw_lo_lo  = nan(1,7); 
    ccw_ccw_lo_lo = nan(1,7);
    
    cw_cw_lo_hi   = nan(1,7);
    cw_ccw_lo_hi  = nan(1,7);    
    ccw_cw_lo_hi  = nan(1,7);
    ccw_ccw_lo_hi = nan(1,7);
    
    % ==> counts
    % ==> high contrast stimuli (==> <choice>_<context>_<contrast>_<dynamic range>)
    cw_cw_hi_lo   = nan(1,7);
    cw_ccw_hi_lo  = nan(1,7);    
    ccw_cw_hi_lo  = nan(1,7);
    ccw_ccw_hi_lo = nan(1,7);
    
    cw_cw_hi_hi   = nan(1,7);
    cw_ccw_hi_hi  = nan(1,7);    
    ccw_cw_hi_hi  = nan(1,7);
    ccw_ccw_hi_hi = nan(1,7);    
    
    % ==> proportions cw for lo contrast
    prop_cw_lcr_ldyn  = nan(1,7);
    prop_ccw_lcr_ldyn = nan(1,7);
    prop_cw_hcr_ldyn  = nan(1,7);
    prop_ccw_hcr_ldyn = nan(1,7);
    
    prop_cw_lcr_hdyn  = nan(1,7);
    prop_ccw_lcr_hdyn = nan(1,7);
    prop_cw_hcr_hdyn  = nan(1,7);
    prop_ccw_hcr_hdyn = nan(1,7);   
    
    % ==> compute counts for each stimulus orientation
    for th = 1:length(or)
        
        % ==> low contrast trials only ==================================================================================================================
        % ==> low dynamic range split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
        if strcmp(dyn_type,'stim')
            % ==> dynamic range for session % => indexing variables - dynamic range
            dv = dvs((ctx == 1) & (ctr == min(ctr)) & ori==or(th),:);
            dynr = dynf(dv); ldyn=(dynf(dvs) <  median(dynr)); hdyn=(dynf(dvs) >= median(dynr));              DV_dynr_meds_check(iS,th,2,1) = median(dynr);    
        end
        % ==> cw congruent choice, lo ctr, lo dynr
        cw_cw_lo_lo(th)   = sum(cw_lcr_cng   & ori==or(th) & ldyn);
        % ==> cw congruent choice, lo ctr, hi dynr
        cw_cw_lo_hi(th)   = sum(cw_lcr_cng   & ori==or(th) & hdyn);                
        % ==> ccw incongruent choice, lo ctr, lo dynr
        ccw_cw_lo_lo(th)  = sum(cw_lcr_icg   & ori==or(th) & ldyn);
        % ==> ccw incongruent choice, lo ctr, hi dynr
        ccw_cw_lo_hi(th)  = sum(cw_lcr_icg   & ori==or(th) & hdyn);                        
               
        if strcmp(dyn_type,'stim')        
            % ==> dynamic range for session % => indexing variables - dynamic range
            dv = dvs((ctx == -1) & (ctr == min(ctr)) & ori==or(th),:);
            dynr = dynf(dv); ldyn=(dynf(dvs) <  median(dynr)); hdyn=(dynf(dvs) >= median(dynr));              DV_dynr_meds_check(iS,th,1,1) = median(dynr);  
        end
        % ==> cw incongruent choice, lo ctr, lo dynr
        cw_ccw_lo_lo(th)  = sum(ccw_lcr_icg  & ori==or(th) & ldyn);    
        % ==> cw incongruent choice, lo ctr, hi dynr
        cw_ccw_lo_hi(th)  = sum(ccw_lcr_icg  & ori==or(th) & hdyn);                                  
        % ==> ccw congruent choice, lo ctr, lo dynr
        ccw_ccw_lo_lo(th) = sum(ccw_lcr_cng  & ori==or(th) & ldyn);      
        % ==> ccw congruent choice, lo ctr, hi dynr
        ccw_ccw_lo_hi(th) = sum(ccw_lcr_cng  & ori==or(th) & hdyn);        

    
        % ==> high contrast trials only =================================================================================================================
        % ==> low dynamic range split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if strcmp(dyn_type,'stim') 
            % ==> dynamic range for session % => indexing variables - dynamic range
            dv = dvs((ctx == 1) & (ctr == max(ctr)) & ori==or(th),:);
            dynr = dynf(dv); ldyn=(dynf(dvs) <  median(dynr)); hdyn=(dynf(dvs) >= median(dynr));              DV_dynr_meds_check(iS,th,2,2) = median(dynr);            
        end
        % ==> cw congruent choice, hi ctr, lo dynr
        cw_cw_hi_lo(th)   = sum(cw_hcr_cng   & ori==or(th) & ldyn);
        % ==> cw congruent choice, hi ctr, lo dynr
        cw_cw_hi_hi(th)   = sum(cw_hcr_cng   & ori==or(th) & hdyn);                                       
        % ==> ccw incongruent choice, hi ctr, lo dynr
        ccw_cw_hi_lo(th)  = sum(cw_hcr_icg   & ori==or(th) & ldyn);
        % ==> ccw incongruent choice, hi ctr, lo dynr
        ccw_cw_hi_hi(th)  = sum(cw_hcr_icg   & ori==or(th) & hdyn);   
            
        if strcmp(dyn_type,'stim')         
            % ==> dynamic range for session % => indexing variables - dynamic range
            dv = dvs((ctx == -1) & (ctr == max(ctr)) & ori==or(th),:);
            dynr = dynf(dv); ldyn=(dynf(dvs) <  median(dynr)); hdyn=(dynf(dvs) >= median(dynr));              DV_dynr_meds_check(iS,th,1,2) = median(dynr);            
        end
        % ==> cw incongruent choice, hi ctr, lo dynr
        cw_ccw_hi_lo(th)  = sum(ccw_hcr_icg  & ori==or(th) & ldyn); 
        % ==> cw incongruent choice, hi ctr, lo dynr
        cw_ccw_hi_hi(th)  = sum(ccw_hcr_icg  & ori==or(th) & hdyn);                       
        % ==> ccw congruent choice, hi ctr, lo dynr
        ccw_ccw_hi_lo(th) = sum(ccw_hcr_cng  & ori==or(th) & ldyn);
        % ==> ccw congruent choice, hi ctr, lo dynr
        ccw_ccw_hi_hi(th) = sum(ccw_hcr_cng  & ori==or(th) & hdyn);                  

        
        % ==> compute ratios for actual data
        prop_cw_fun();      
    end
    
    % ==> PF curve fit
    fprintf('computing psychometric functions...\n')    
    % ==> low dynamic range cases (low contrast stimulus trials)
    [PF_ldyn_lcr, p_ldyn_lcr] = PF_fit_fun(ccw_cw_lo_lo, ccw_ccw_lo_lo, cw_cw_lo_lo, cw_ccw_lo_lo, or, startVec_M1, LB_M1, UB_M1, options);
    % ==> high dynamic range cases (low contrast stimulus trials)
    [PF_hdyn_lcr, p_hdyn_lcr] = PF_fit_fun(ccw_cw_lo_hi, ccw_ccw_lo_hi, cw_cw_lo_hi, cw_ccw_lo_hi, or, startVec_M1, LB_M1, UB_M1, options);
    
    % ==> low dynamic range cases (high contrast stimulus trials)
    [PF_ldyn_hcr, p_ldyn_hcr] = PF_fit_fun(ccw_cw_hi_lo, ccw_ccw_hi_lo, cw_cw_hi_lo, cw_ccw_hi_lo, or, startVec_M1, LB_M1, UB_M1, options);
    % ==> high dynamic range cases (high contrast stimulus trials)        
    [PF_hdyn_hcr, p_hdyn_hcr] = PF_fit_fun(ccw_cw_hi_hi, ccw_ccw_hi_hi, cw_cw_hi_hi, cw_ccw_hi_hi, or, startVec_M1, LB_M1, UB_M1, options);
    
    % ==> plots
    figure(3);    
    ax = subplot(5,6,iS); 
    title([S.general.monkey, S.general.expDate]); %title(['Session: ', num2str(iS)]);
    hold on; hold all;
    plot(or,PF_ldyn_lcr(1,:),'b.-','linewidth', 1)    
    plot(or,PF_ldyn_lcr(2,:),'r.-','linewidth', 1)    
    plot(or,PF_hdyn_lcr(1,:),'b:','linewidth', 1)
    plot(or,PF_hdyn_lcr(2,:),'r:','linewidth', 1)   
    plot(zeros(length(or),1),linspace(0,1,length(or)),'k--','linewidth', 0.5) 
    plot(or,ones(length(or),1)/2,'k--','linewidth', 0.5) 
    ax.XTick = or;
    ax.XTickLabel = or;
    xlabel('orientation')
    ylabel('prop cw')
    drawnow;
    
    % ==> low dynamic range (and low contrast)
    % ==> count trials in ccw condition (==> <choice>_<context>_<contrast>_<dynamic range>)
    tn_ccw_lo = cw_ccw_lo_lo + ccw_ccw_lo_lo;
    % ==> count trials in cw condition
    tn_cw_lo = cw_cw_lo_lo + ccw_cw_lo_lo;    
    % ==> high dynamic range (and low contrast)
    % ==> count trials in ccw condition (==> <choice>_<context>_<contrast>_<dynamic range>)
    tn_ccw_hi = cw_ccw_lo_hi + ccw_ccw_lo_hi;
    % ==> count trials in cw condition
    tn_cw_hi = cw_cw_lo_hi + ccw_cw_lo_hi;    
    
    % ==> trial proportions for plot (just count, for each stimulus type, number of trials over total session trials)
    tn_ccw_lo_p = tn_ccw_lo ./ sum(tn_ccw_lo + tn_cw_lo); %ccw_ccw_lo_lo ./ sum(tn_ccw_lo);
    tn_cw_lo_p =  tn_cw_lo  ./ sum(tn_ccw_lo + tn_cw_lo); %cw_cw_lo_lo   ./ sum(tn_cw_lo); 
    % ==> high dynamic range
    tn_ccw_hi_p = tn_ccw_hi ./ sum(tn_ccw_hi + tn_cw_hi); %ccw_ccw_lo_hi ./ sum(tn_ccw_hi);
    tn_cw_hi_p =  tn_cw_hi  ./ sum(tn_ccw_hi + tn_cw_hi);  %cw_cw_lo_hi   ./ sum(tn_cw_hi);  
    
    % ==> plot trial proportions
    for o = 1:length(or)
        hold on; hold all;
        % ==> trial proportions for low dynamic range
        scatter(or(o), PF_ldyn_lcr(2,o),  max([1,round(tn_ccw_lo_p(o)*1000)]), 'o','filled', 'markerfacecolor', [1,0,0], 'markeredgecolor','r')
        scatter(or(o), PF_ldyn_lcr(1,o),  max([1,round(tn_cw_lo_p(o)*1000)]),  'o','filled', 'markerfacecolor', [0,0,1], 'markeredgecolor','b')
        % ==> trial proportions for high dynamic range
        scatter(or(o), PF_hdyn_lcr(2,o),  max([1,round(tn_ccw_hi_p(o)*1000)]), 'o','filled', 'markerfacecolor', [1,0.5,0.5], 'markeredgecolor','r')
        scatter(or(o), PF_hdyn_lcr(1,o),  max([1,round(tn_cw_hi_p(o)*1000)]),  'o','filled', 'markerfacecolor', [0.5,0.5,1], 'markeredgecolor','b')
    end  
    axis square;
    drawnow;    

%     % ==> plots
%     figure(4);
%     subplot(5,6,iS); title(['Session: ', num2str(iS)]);
%     plot(or,prop_cw_lcr_hdyn,'b.-','linewidth', 1)
%     hold on; hold all;
%     plot(or,prop_ccw_lcr_hdyn,'r.-','linewidth', 1)
%     hold on; hold all;
%     plot(or,prop_cw_lcr_ldyn,'b+:','linewidth', 1)
%     hold on; hold all;
%     plot(or,prop_ccw_lcr_ldyn,'rx:','linewidth', 1)   
%     xlabel('orientation')
%     ylabel('prop cw')
%     drawnow;    
%     
%     % ==> plots
%     figure(5);
%     subplot(5,6,iS); title(['Session: ', num2str(iS)]);
%     plot(or,PF_ldyn_hcr(1,:),'g.-','linewidth', 1)
%     hold on; hold all;
%     plot(or,PF_ldyn_hcr(2,:),'m.-','linewidth', 1)
%     hold on; hold all;
%     plot(or,PF_hdyn_hcr(1,:),'g:','linewidth', 1)
%     hold on; hold all;
%     plot(or,PF_hdyn_hcr(2,:),'m:','linewidth', 1)   
%     xlabel('orientation')
%     ylabel('prop cw')
%     drawnow;
%     
%     figure(6);
%     subplot(5,6,iS); title(['Session: ', num2str(iS)]);    
%     plot(or,prop_cw_hcr_hdyn,'g+:')
%     hold on; hold all;
%     plot(or,prop_ccw_hcr_hdyn,'mx:')
%     hold on; hold all;       
%     plot(or,prop_cw_hcr_ldyn,'g.-')
%     hold on; hold all;
%     plot(or,prop_ccw_hcr_ldyn,'m.-')
%     hold on; hold all;    
%     xlabel('orientation')
%     ylabel('prop cw')
%     drawnow;    
    
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

% ==> quick check
if strcmp(dyn_type,'stim') 
    assert(unique(DV_dynr_meds_check == DV_dynr_meds))
end

% ==> read in important metrics for final analysis
dvCatPerf = load([drc,'dvCatPerf.mat']);
dvCatPerf = dvCatPerf.dvCatPerf';

% ==> AIC Delta (from DV peak glm analysis in F2)
aic_d = load([drc,'aic_d.mat']);
aic_d = aic_d.aic_d';

% ==> log-likelihood ratios (from DV peak glm analysis in F2)
llr = load([drc,'llr.mat']);
llr = llr.llr';

%%
close all; clc;
ppl = [1:13]; ids = [1:29];
fg = figure(13); clf; set(fg,'color','white'); set(fg,'Position',[723 4 1104 958]);
pltf([db_lcr], [dp_lcr], [aic_d], [llr], [dvCatPerf], ppl, 'lo_hi_ctr', ids);

ppl = [1:13]; ids = [1:29];
fg = figure(14); clf; set(fg,'color','white'); set(fg,'Position',[723 4 1104 958]);
pltf([db_hcr], [dp_hcr], [aic_d], [llr], [dvCatPerf], ppl, 'lo_hi_ctr', ids);

ppl = [14:29]; ids = [1:29];
fg = figure(15); clf; set(fg,'color','white'); set(fg,'Position',[723 4 1104 958]);
pltf([db_lcr], [dp_lcr], [aic_d], [llr], [dvCatPerf], ppl, 'lo_hi_ctr', ids);

ppl = [14:29]; ids = [1:29];
fg = figure(16); clf; set(fg,'color','white'); set(fg,'Position',[723 4 1104 958]);
pltf([db_hcr], [dp_hcr], [aic_d], [llr], [dvCatPerf], ppl, 'lo_hi_ctr', ids);

%% ==> DV choice predictivity, AIC delta analysis etc

% close all; clc;
ppl = [1:29,(1:29)+29]; ids = [1:29,1:29];
fg = figure(12); clf; set(fg,'color','white'); set(fg,'Position',[723 4 1104 958]);
pltf([db_lcr;db_hcr], [dp_lcr;dp_hcr], [aic_d;aic_d], [llr;llr], [dvCatPerf;dvCatPerf], ppl, 'lo_hi_ctr', ids);

%%
% ==> F
% close all; clc;
ppl = [1:13,(1:13)+29]; ids = [1:29,1:29];
fg = figure(13); clf; set(fg,'color','white'); set(fg,'Position',[723 4 1104 958]);
pltf([db_lcr;db_hcr], [dp_lcr;dp_hcr], [aic_d;aic_d], [llr;llr], [dvCatPerf;dvCatPerf], ppl, 'lo_hi_ctr', ids);

% ==> JP
% close all; clc;
ppl = [14:29,(14:29)+29]; ids = [1:29,1:29];
fg = figure(15); clf; set(fg,'color','white'); set(fg,'Position',[723 4 1104 958]);
pltf([db_lcr;db_hcr], [dp_lcr;dp_hcr], [aic_d;aic_d], [llr;llr], [dvCatPerf;dvCatPerf], ppl, 'lo_hi_ctr', ids);

%% ==> Average delta bias by animal and by contrast (F4G)
clc;

% ==> plot average delta Bias
figure(17); set(gcf,'color','white');
hold on; hold all;
% ==> monkey F - lo & hi contrast
errorbar([1,2], ...
         [mean(db_lcr(1:13)), mean(db_hcr(1:13))], ...
         [std(db_lcr(1:13))/sqrt(13), std(db_hcr(1:13))/sqrt(13)] , ...
         'o','color','k');
% ==> monkey J - lo & hi contrast
errorbar([3,4], ...
         [mean(db_lcr(14:end)), mean(db_hcr(14:end))], ...
         [std(db_lcr(14:end))/sqrt(16), std(db_hcr(14:end))/sqrt(16)] , ...
         'ko');     
% ==> lo & hi contrast     
errorbar([5,6], ...
         [mean(db_lcr), mean(db_hcr)], ...
         [std(db_lcr)/sqrt(length(db_lcr)), std(db_hcr)/sqrt(length(db_hcr))] , ...
         'ko');
% ==> all
errorbar([7], ...
         [mean([db_lcr;db_hcr])], ...
         [std([db_lcr;db_hcr])/sqrt(length(db_lcr)*2)], ...
         'ko');
xlim([0,7]);


% ==> wilcoxon tests
% => low contrast cases
[pJl,~] = signrank(db_lcr(14:end))
[pFl,~] = signrank(db_lcr(1:13))

% => hicg contrast cases
[pJh,~] = signrank(db_hcr(14:end))
[pFh,~] = signrank(db_hcr(1:13))

[pl,~] = signrank(db_lcr)
[ph,~] = signrank(db_hcr)

[pall,~] = signrank([db_lcr; db_hcr])

%% ==> correlation plots (delta bias and delta AIC) (F4H)
clc;

rFl = corr(db_lcr(1:13),aic_d(1:13));
rFh = corr(db_hcr(1:13),aic_d(1:13));

rJl = corr(db_lcr(14:end),aic_d(14:end));
rJh = corr(db_hcr(14:end),aic_d(14:end));

% ==> fisher transformation function
z = @(r) (1/2)*log((1 + r)/(1 - r));
% ==> inverse of the z transformation
zi = @(z) (exp(2*z) - 1)/(exp(2*z) + 1);

zFli = zi(z(rFl));
% ==> should match
% assert(zFli == tanh(atanh(rFl)));
zFhi = zi(z(rFh));
% ==> should match
% assert(zFhi == tanh(atanh(rFh)));

zJli = zi(z(rJl));
% ==> should match
% assert(zJli == tanh(atanh(rJl)));
zJhi = zi(z(rJh));
% ==> should match
% assert(zJhi == tanh(atanh(rJh)));

% => standard deviation (1/sqrt(N - 3))
% => monkey F (N = 13), so SD denom = sqrt(N-3)
sdFl_ub = z(rFl) + 1/sqrt(10); sdFl_lb = z(rFl) - 1/sqrt(10); 
sdFh_ub = z(rFh) + 1/sqrt(10); sdFh_lb = z(rFh) - 1/sqrt(10);
% => monkey J (N = 16), so SD denom = sqrt(N-3)
sdJl_ub = z(rJl) + 1/sqrt(13); sdJl_lb = z(rJl) - 1/sqrt(13); 
sdJh_ub = z(rJh) + 1/sqrt(13); sdJh_lb = z(rJh) - 1/sqrt(13);

%NOTE: (errorbar(X,Y,NEG(i),POS(i)))
% ==> errorbar plot allows NEG(i), POS(i), which are ~distances~ above and below
% correlations for drawing error bars in the plot. So these are computed as follows:

% => NEG(i) is correlation r - z^(-1)(SD_LB), where SD_LB is the lower
% bound correlation in z coordinates, e.g. 
% ==> SD_LB = z(r) - 1/sqrt(N-3)
% so the ~distance~ below r for NEG(i) is r - z^(-1)(SD_LB)
% (actual correlation minus the lower bound correlation)

% => POS(i) is correlation z(-1)(SD_UB), where SD_UB is the higher 
% bound correlation in z coordinates, e.g.
% ==> SD_UB = r(z) + 1/sqrt(N-3)
% so the ~distance~ above r for POS(i) is z(-1)(SD_UB) - r
% (upper bound correlation minus actual correlation)

figure(18); set(gcf,'color','white');
hold on; hold all;
% ==> monkey F 
errorbar([1,2], ...
         [rFl, rFh], ... % ==> original correlations
         [rFl-zi(sdFl_lb), rFh-zi(sdFh_lb)],[zi(sdFl_ub)-rFl, zi(sdFh_ub)-rFh] , ...
         'o','color','k');
% ==> monkey J
errorbar([3,4], ...
         [rJl, rJh], ... % ==> original correlations
         [rJl-zi(sdJl_lb), rJh-zi(sdJh_lb)], [zi(sdJl_ub)-rJl, zi(sdJh_ub)-rJh] , ...
         'ko'); 
xlim([0,7]);
ylim([0,1]);

% ==> z-scores across animals

db_lcr_z = nan(1,29);
db_hcr_z = nan(1,29);

% ==> Animal F
% ==> low contrast
db_lcr_z(1:13) = (db_lcr(1:13) - mean(db_lcr(1:13)))/std(db_lcr(1:13));
% ==> high contrast
db_hcr_z(1:13) = (db_hcr(1:13) - mean(db_hcr(1:13)))/std(db_hcr(1:13)); 
% ==> Animal J
% ==> low contrast
db_lcr_z(14:29) = (db_lcr(14:29) - mean(db_lcr(14:29)))/std(db_lcr(14:29));
% ==> high contrast
db_hcr_z(14:29) = (db_hcr(14:29) - mean(db_hcr(14:29)))/std(db_hcr(14:29)); 

% ==> standardize aic by monkey

aic_d_z(1:13)  = (aic_d(1:13)  - mean(aic_d(1:13)))/std(aic_d(1:13));
aic_d_z(14:29) = (aic_d(14:29) - mean(aic_d(14:29)))/std(aic_d(14:29));

% ==> correlation
[r_z_lcr, p_z_lcr] = corr(aic_d_z',db_lcr_z');
[r_z_hcr, p_z_hcr] = corr(aic_d_z',db_hcr_z');

% ==> correlation across all contrasts and animals
[rz_db, pz_db] = corr([aic_d_z,aic_d_z]',[db_lcr_z, db_hcr_z]');

lcr_se_lb = z(r_z_lcr) - 1/sqrt(26); 
lcr_se_ub = z(r_z_lcr) + 1/sqrt(26);

hcr_se_lb = z(r_z_hcr) - 1/sqrt(26); 
hcr_se_ub = z(r_z_hcr) + 1/sqrt(26);

se_lb = z(rz_db) - 1/sqrt(29*2 - 3); 
se_ub = z(rz_db) + 1/sqrt(29*2 - 3);

figure(18);
errorbar([5,6,7], ...
         [r_z_lcr, r_z_hcr, rz_db], ...
         [r_z_lcr-zi(lcr_se_lb), r_z_hcr-zi(hcr_se_lb), rz_db-zi(se_lb)], ...
         [zi(lcr_se_ub)-r_z_lcr, zi(hcr_se_ub)-r_z_hcr, zi(se_ub)-rz_db], ...
         'ro');
ylabel('\Delta bias and \Delta AIC')
     
figure(19); 
subplot(1,2,1);
set(gcf,'color','white');
scatter(aic_d_z,db_lcr_z,50, 'ko', 'filled');
title(['r = ', num2str(r_z_lcr), ', p = ',num2str(p_z_lcr)]);
ylabel('Standardized \Delta bias (z)');
xlabel('Standardized \Delta AIC (z)');
axis square;
subplot(1,2,2);
scatter(aic_d_z,db_hcr_z,50, 'ko', 'filled');
title(['r = ', num2str(r_z_hcr), ', p = ',num2str(p_z_hcr)]);
ylabel('Standardized \Delta bias (z)');
xlabel('Standardized \Delta AIC (z)');
axis square;     

%% ==> average delta perceptual uncertainty (F4I)
clc;

% ==> plot average delta perceptual uncertainty
figure(20); set(gcf,'color','white');
hold on; hold all;
% ==> monkey F
errorbar([1,2], ...
         [mean(dp_lcr(1:13)), mean(dp_hcr(1:13))], ...
         [std(dp_lcr(1:13))/sqrt(13), std(dp_hcr(1:13))/sqrt(13)] , ...
         'o','color','k');
% ==> monkey J
errorbar([3,4], ...
         [mean(dp_lcr(14:end)), mean(dp_hcr(14:end))], ...
         [std(dp_lcr(14:end))/sqrt(16), std(dp_hcr(14:end))/sqrt(16)] , ...
         'ko');     
% ==> lo & hi contrast     
errorbar([5,6], ...
         [mean(dp_lcr), mean(dp_hcr)], ...
         [std(dp_lcr)/sqrt(length(dp_lcr)), std(dp_hcr)/sqrt(length(dp_hcr))] , ...
         'ko');
% ==> all
errorbar([7], ...
         [mean([dp_lcr;dp_hcr])], ...
         [std([dp_lcr;dp_hcr])/sqrt(length(dp_lcr)*2)], ...
         'ko');
xlim([0,7]);     


% ==> wilcoxon tests
% => low contrast cases
[pJlu,~] = signrank(dp_lcr(14:end))
[pFlu,~] = signrank(dp_lcr(1:13))

% => hicg contrast cases
[pJhu,~] = signrank(dp_hcr(14:end))
[pFhu,~] = signrank(dp_hcr(1:13))

[pl,~] = signrank(dp_lcr)
[ph,~] = signrank(dp_hcr)

[pall,~] = signrank([dp_lcr; dp_hcr])

%% ==> correlation delta perceptual uncertainty and delta AIC (F4J)

% ==> make a function from this spaghetti

rFlpdpu = corr(dp_lcr(1:13),aic_d(1:13));
rFhpdpu = corr(dp_hcr(1:13),aic_d(1:13));

rJlpdpu = corr(dp_lcr(14:end),aic_d(14:end));
rJhpdpu = corr(dp_hcr(14:end),aic_d(14:end));

zFli = zi(z(rFlpdpu));
% ==> should match
% assert(zFli == tanh(atanh(rFl)));
zFhi = zi(z(rFhpdpu));
% ==> should match
% assert(zFhi == tanh(atanh(rFh)));

zJli = zi(z(rJlpdpu));
% ==> should match
% assert(zJli == tanh(atanh(rJl)));
zJhi = zi(z(rJhpdpu));
% ==> should match
% assert(zJhi == tanh(atanh(rJh)));

% => standard deviation (1/sqrt(N - 3))
% => monkey F (N = 13), so SD denom = sqrt(N-3)
sdFl_ub = z(rFlpdpu) + 1/sqrt(10); sdFl_lb = z(rFlpdpu) - 1/sqrt(10); 
sdFh_ub = z(rFhpdpu) + 1/sqrt(10); sdFh_lb = z(rFhpdpu) - 1/sqrt(10);
% => monkey J (N = 16), so SD denom = sqrt(N-3)
sdJl_ub = z(rJlpdpu) + 1/sqrt(13); sdJl_lb = z(rJlpdpu) - 1/sqrt(13); 
sdJh_ub = z(rJhpdpu) + 1/sqrt(13); sdJh_lb = z(rJhpdpu) - 1/sqrt(13);

%NOTE: (errorbar(X,Y,NEG(i),POS(i)))
% ==> errorbar plot allows NEG(i), POS(i), which are ~distances~ above and below
% correlations for drawing error bars in the plot. So these are computed as follows:

% => NEG(i) is correlation r - z^(-1)(SD_LB), where SD_LB is the lower
% bound correlation in z coordinates, e.g. 
% ==> SD_LB = z(r) - 1/sqrt(N-3)
% so the ~distance~ below r for NEG(i) is r - z^(-1)(SD_LB)
% (actual correlation minus the lower bound correlation)

% => POS(i) is correlation z(-1)(SD_UB), where SD_UB is the higher 
% bound correlation in z coordinates, e.g.
% ==> SD_UB = r(z) + 1/sqrt(N-3)
% so the ~distance~ above r for POS(i) is z(-1)(SD_UB) - r
% (upper bound correlation minus actual correlation)

figure(21); set(gcf,'color','white');
hold on; hold all;
% ==> monkey F 
errorbar([1,2], ...
         [rFlpdpu, rFhpdpu], ... % ==> original correlations
         [rFlpdpu-zi(sdFl_lb), rFhpdpu-zi(sdFh_lb)],[zi(sdFl_ub)-rFlpdpu, zi(sdFh_ub)-rFhpdpu] , ...
         'o','color','k');
% ==> monkey J
errorbar([3,4], ...
         [rJlpdpu, rJhpdpu], ... % ==> original correlations
         [rJlpdpu-zi(sdJl_lb), rJhpdpu-zi(sdJh_lb)], [zi(sdJl_ub)-rJlpdpu, zi(sdJh_ub)-rJhpdpu] , ...
         'ko'); 
xlim([0,7]);
ylim([0,1]);

% ==> z-scores across animals

dp_lcr_z = nan(1,29);
dp_hcr_z = nan(1,29);

% ==> Animal F
% ==> low contrast
dp_lcr_z(1:13) = (dp_lcr(1:13) - mean(dp_lcr(1:13)))/std(dp_lcr(1:13));
% ==> high contrast
dp_hcr_z(1:13) = (dp_hcr(1:13) - mean(dp_hcr(1:13)))/std(dp_hcr(1:13)); 
% ==> Animal J
% ==> low contrast
dp_lcr_z(14:29) = (dp_lcr(14:29) - mean(dp_lcr(14:29)))/std(dp_lcr(14:29));
% ==> high contrast
dp_hcr_z(14:29) = (dp_hcr(14:29) - mean(dp_hcr(14:29)))/std(dp_hcr(14:29)); 

% ==> standardize aic by monkey

aic_d_z(1:13)  = (aic_d(1:13)  - mean(aic_d(1:13)))/std(aic_d(1:13));
aic_d_z(14:29) = (aic_d(14:29) - mean(aic_d(14:29)))/std(aic_d(14:29));

% ==> correlation
[r_z_lcr, p_z_lcr] = corr(aic_d_z',dp_lcr_z');
[r_z_hcr, p_z_hcr] = corr(aic_d_z',dp_hcr_z');

% ==> correlation across all contrasts and animals
[rz_dp, pz_dp] = corr([aic_d_z,aic_d_z]',[dp_lcr_z, dp_hcr_z]');

lcr_se_lb = z(r_z_lcr) - 1/sqrt(26); 
lcr_se_ub = z(r_z_lcr) + 1/sqrt(26);

hcr_se_lb = z(r_z_hcr) - 1/sqrt(26); 
hcr_se_ub = z(r_z_hcr) + 1/sqrt(26);

se_lb = z(rz_dp) - 1/sqrt(29*2 - 3); 
se_ub = z(rz_dp) + 1/sqrt(29*2 - 3);

errorbar([5,6,7], ...
         [r_z_lcr, r_z_hcr, rz_dp], ...
         [r_z_lcr-zi(lcr_se_lb), r_z_hcr-zi(hcr_se_lb), rz_dp-zi(se_lb)], ... 
         [zi(lcr_se_ub)-r_z_lcr, zi(hcr_se_ub)-r_z_hcr, zi(se_ub)-rz_dp], ...
         'ro');
     
figure(22); 
subplot(1,2,1);
set(gcf,'color','white');
scatter(aic_d_z,dp_lcr_z,50, 'ko', 'filled');
title(['r = ', num2str(r_z_lcr), ', p = ',num2str(p_z_lcr)]);
ylabel('Standardized \Delta perceptual uncertainty (z)');
xlabel('Standardized \Delta AIC (z)');
axis square;
subplot(1,2,2);
scatter(aic_d_z,dp_hcr_z,50, 'ko', 'filled');
title(['r = ', num2str(r_z_hcr), ', p = ',num2str(p_z_hcr)]);
ylabel('Standardized \Delta perceptual uncertainty (z)');
xlabel('Standardized \Delta AIC (z)');
axis square;

%% ==> correlations between db and dp

% => create a function for this spaghetti (function just makes plots with z and z^(-1) transforms)

rFlpdpu = corr(db_lcr(1:13), dp_lcr(1:13));
rFhpdpu = corr(db_hcr(1:13), dp_hcr(1:13));

rJlpdpu = corr(db_lcr(14:end), dp_lcr(14:end));
rJhpdpu = corr(db_hcr(14:end), dp_hcr(14:end));

zFli = zi(z(rFlpdpu));
zFhi = zi(z(rFhpdpu));

zJli = zi(z(rJlpdpu));
zJhi = zi(z(rJhpdpu));

% => standard deviation (1/sqrt(N - 3))
% => monkey F (N = 13), so SD denom = sqrt(N-3)
sdFl_ub = z(rFlpdpu) + 1/sqrt(10); sdFl_lb = z(rFlpdpu) - 1/sqrt(10); 
sdFh_ub = z(rFhpdpu) + 1/sqrt(10); sdFh_lb = z(rFhpdpu) - 1/sqrt(10);
% => monkey J (N = 16), so SD denom = sqrt(N-3)
sdJl_ub = z(rJlpdpu) + 1/sqrt(13); sdJl_lb = z(rJlpdpu) - 1/sqrt(13); 
sdJh_ub = z(rJhpdpu) + 1/sqrt(13); sdJh_lb = z(rJhpdpu) - 1/sqrt(13);

%NOTE: (errorbar(X,Y,NEG(i),POS(i)))
% ==> errorbar plot allows NEG(i), POS(i), which are ~distances~ above and below
% correlations for drawing error bars in the plot. So these are computed as follows:

% => NEG(i) is correlation r - z^(-1)(SD_LB), where SD_LB is the lower
% bound correlation in z coordinates, e.g. 
% ==> SD_LB = z(r) - 1/sqrt(N-3)
% so the ~distance~ below r for NEG(i) is r - z^(-1)(SD_LB)
% (actual correlation minus the lower bound correlation)

% => POS(i) is correlation z(-1)(SD_UB), where SD_UB is the higher 
% bound correlation in z coordinates, e.g.
% ==> SD_UB = r(z) + 1/sqrt(N-3)
% so the ~distance~ above r for POS(i) is z(-1)(SD_UB) - r
% (upper bound correlation minus actual correlation)

figure(21); set(gcf,'color','white');
hold on; hold all;
% ==> monkey F 
errorbar([1,2], ...
         [rFlpdpu, rFhpdpu], ... % ==> original correlations
         [rFlpdpu-zi(sdFl_lb), rFhpdpu-zi(sdFh_lb)],[zi(sdFl_ub)-rFlpdpu, zi(sdFh_ub)-rFhpdpu] , ...
         'o','color','k');
% ==> monkey J
errorbar([3,4], ...
         [rJlpdpu, rJhpdpu], ... % ==> original correlations
         [rJlpdpu-zi(sdJl_lb), rJhpdpu-zi(sdJh_lb)], [zi(sdJl_ub)-rJlpdpu, zi(sdJh_ub)-rJhpdpu] , ...
         'ko'); 
xlim([0,7]);
ylim([-0.1,1]);

% ==> z-scores across animals

% ==> correlation
[r_z_lcr, p_z_lcr] = corr(db_lcr_z',dp_lcr_z');
[r_z_hcr, p_z_hcr] = corr(db_hcr_z',dp_hcr_z');

% ==> correlation across all contrasts and animals
[rz_dp, pz_dp] = corr([db_lcr_z,db_hcr_z]',[dp_lcr_z, dp_hcr_z]');

lcr_se_lb = z(r_z_lcr) - 1/sqrt(26); 
lcr_se_ub = z(r_z_lcr) + 1/sqrt(26);

hcr_se_lb = z(r_z_hcr) - 1/sqrt(26); 
hcr_se_ub = z(r_z_hcr) + 1/sqrt(26);

se_lb = z(rz_dp) - 1/sqrt(29*2 - 3); 
se_ub = z(rz_dp) + 1/sqrt(29*2 - 3);

errorbar([5,6,7], ...
         [r_z_lcr, r_z_hcr, rz_dp], ...
         [r_z_lcr-zi(lcr_se_lb), r_z_hcr-zi(hcr_se_lb), rz_dp-zi(se_lb)], ... 
         [zi(lcr_se_ub)-r_z_lcr, zi(hcr_se_ub)-r_z_hcr, zi(se_ub)-rz_dp], ...
         'ro');
     
figure(22); 
subplot(1,2,1);
set(gcf,'color','white');
scatter(aic_d_z,dp_lcr_z,50, 'ko', 'filled');
title(['r = ', num2str(r_z_lcr), ', p = ',num2str(p_z_lcr)]);
ylabel('Standardized \Delta perceptual uncertainty (z)');
xlabel('Standardized \Delta AIC (z)');
axis square;
subplot(1,2,2);
scatter(aic_d_z,dp_hcr_z,50, 'ko', 'filled');
title(['r = ', num2str(r_z_hcr), ', p = ',num2str(p_z_hcr)]);
ylabel('Standardized \Delta perceptual uncertainty (z)');
xlabel('Standardized \Delta AIC (z)');
axis square;

%% ==> example plot

figure(); set(gcf,'color','white');
scatter(db_lcr(14:29), dp_lcr(14:29),75,'ko','filled');
axis square;

%% ==> ANCOVA

dbz  = [db_lcr_z, db_hcr_z];
dpz  = [dp_lcr_z, dp_hcr_z];
aicz = [aic_d_z, aic_d_z];

% ==> median split to return aic grouping variable (gv)
% aicgv = aicz >= median(aicz);

% => histogram to get grouping variable.
[~,~,aicgv] = histcounts(aicz,8);

% ==> run ANCOVA
aoctool(dbz,dpz,aicgv)

