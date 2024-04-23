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
% dyn_type = 'popu';
% ==> estimate dynamic range per stimulus?
dyn_type = 'stim';

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
figure(4); set(gcf,'color','white'); set(gcf,'Position',[86 282 921 668]);
figure(5); set(gcf,'color','white'); set(gcf,'Position',[86 282 921 668]);
figure(6); set(gcf,'color','white'); set(gcf,'Position',[86 282 921 668]);


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
    % ==> lo contrast stimuli
    cw_cw_lo_lo   = nan(1,7); 
    cw_ccw_lo_lo  = nan(1,7);    
    ccw_cw_lo_lo  = nan(1,7); %cho_ccw_ct_cw_cr_lo_dynr_lo
    ccw_ccw_lo_lo = nan(1,7);
    
    cw_cw_lo_hi   = nan(1,7);
    cw_ccw_lo_hi  = nan(1,7);    
    ccw_cw_lo_hi  = nan(1,7);
    ccw_ccw_lo_hi = nan(1,7);
    
    % ==> counts
    % ==> high contrast stimuli
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
    subplot(5,6,iS); title(['Session: ', num2str(iS)]);
    plot(or,PF_ldyn_lcr(1,:),'b.-','linewidth', 1)
    hold on; hold all;
    plot(or,PF_ldyn_lcr(2,:),'r.-','linewidth', 1)
    hold on; hold all;
    plot(or,PF_hdyn_lcr(1,:),'b+:','linewidth', 1)
    hold on; hold all;
    plot(or,PF_hdyn_lcr(2,:),'rx:','linewidth', 1)   
    xlabel('orientation')
    ylabel('prop cw')
    drawnow;

    % ==> plots
    figure(4);
    subplot(5,6,iS); title(['Session: ', num2str(iS)]);
    plot(or,prop_cw_lcr_hdyn,'b.-','linewidth', 1)
    hold on; hold all;
    plot(or,prop_ccw_lcr_hdyn,'r.-','linewidth', 1)
    hold on; hold all;
    plot(or,prop_cw_lcr_ldyn,'b+:','linewidth', 1)
    hold on; hold all;
    plot(or,prop_ccw_lcr_ldyn,'rx:','linewidth', 1)   
    xlabel('orientation')
    ylabel('prop cw')
    drawnow;
    
%     figure(3);
%     subplot(5,6,iS); 
%     for o = 1:length(or)
%         hold on; hold all;
%         scatter(or(o), prop_cw_lcr_hdyn(o),  max([1,round(prop_cw_lcr_hdyn(o)*100)]), 'o','filled','markerfacecolor', [0,0,1],     'markeredgecolor','b')
%         scatter(or(o), prop_ccw_lcr_hdyn(o), max([1,round(prop_ccw_lcr_hdyn(o)*100)]), 'o','filled','markerfacecolor',[1,0,0],     'markeredgecolor','r')
%         scatter(or(o), prop_cw_lcr_ldyn(o),  max([1,round(prop_cw_lcr_ldyn(o)*100)]), 'o','filled','markerfacecolor', [0.5,0.5,1], 'markeredgecolor','b')
%         scatter(or(o), prop_ccw_lcr_ldyn(o), max([1,round(prop_ccw_lcr_ldyn(o)*100)]),'o','filled','markerfacecolor', [1,0.5,0.5], 'markeredgecolor','r')
%     end
%     drawnow;
    
    % ==> plots
    figure(5);
    subplot(5,6,iS); title(['Session: ', num2str(iS)]);
    plot(or,PF_ldyn_hcr(1,:),'g.-')
    hold on; hold all;
    plot(or,PF_ldyn_hcr(2,:),'m.-')
    hold on; hold all;
    plot(or,PF_hdyn_hcr(1,:),'gx:')
    hold on; hold all;
    plot(or,PF_hdyn_hcr(2,:),'m+:')    
    xlabel('orientation')
    ylabel('prop cw')
    drawnow;
    
    figure(6);
    subplot(5,6,iS); title(['Session: ', num2str(iS)]);    
    plot(or,prop_cw_hcr_hdyn,'g+:')
    hold on; hold all;
    plot(or,prop_ccw_hcr_hdyn,'mx:')
    hold on; hold all;       
    plot(or,prop_cw_hcr_ldyn,'g.-')
    hold on; hold all;
    plot(or,prop_ccw_hcr_ldyn,'m.-')
    hold on; hold all;    
    xlabel('orientation')
    ylabel('prop cw')
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

% ==> plot average delta Bias
figure(17); set(gcf,'color','white');
hold on; hold all;
% ==> monkey F
errorbar([1,2], ...
         [mean(db_lcr(1:13)), mean(db_hcr(1:13))], ...
         [std(db_lcr(1:13))/sqrt(13), std(db_hcr(1:13))/sqrt(13)] , ...
         'o','color','k');
% ==> monkey J
errorbar([3,4], ...
         [mean(db_lcr(14:end)), mean(db_hcr(14:end))], ...
         [std(db_lcr(14:end))/sqrt(16), std(db_hcr(14:end))/sqrt(16)] , ...
         'ko');     
     
xlim([0,5]);

% ==> wilcoxon tests
% => low contrast cases
[pJl,~] = signrank(db_lcr(14:end))
[pFl,~] = signrank(db_lcr(1:13))

% => hich contrast cases
[pJh,~] = signrank(db_hcr(14:end))
[pFh,~] = signrank(db_hcr(1:13))

%% ==> correlation plots (delta bias and delta AIC) (F4H)

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
% => monkey F (N = 13)
sdFl = 1/sqrt(13 - 3); sdFh = 1/sqrt(13 - 3);
% => monkey J (N = 16)
sdJl = 1/sqrt(16 - 3); sdJh = 1/sqrt(16 - 3);

figure(18); set(gcf,'color','white');
hold on; hold all;
% ==> monkey F
errorbar([1,2], ...
         [zFli, zFhi], ...
         [zi(sdFl),zi(sdFh)] , ...
         'o','color','k');
% ==> monkey J
errorbar([3,4], ...
         [zJli, zJhi], ...
         [zi(sdJl),zi(sdJh)] , ...
         'ko');  

% figure(18); set(gcf,'color','white');
% hold on; hold all;
% % ==> monkey F
% errorbar([1,2], ...
%          [z(rFl), z(rFh)], ...
%          [sdFl,sdFh] , ...
%          'o','color','k');
% % ==> monkey J
% errorbar([3,4], ...
%          [z(rJl), z(rJh)], ...
%          [sdJl,sdJh] , ...
%          'ko');     
     
xlim([0,5]);
ylim([0,1.1]);

%% ==> average delta perceptual uncertainty (F4I)

% ==> plot average delta perceptual uncertainty
figure(19); set(gcf,'color','white');
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
     
xlim([0,5]);

% ==> wilcoxon tests
% => low contrast cases
[pJlu,~] = signrank(dp_lcr(14:end))
[pFlu,~] = signrank(dp_lcr(1:13))

% => hich contrast cases
[pJhu,~] = signrank(dp_hcr(14:end))
[pFhu,~] = signrank(dp_hcr(1:13))


%% ==> correlation delta perceptual uncertainty and delta AIC (F4J)

rFl = corr(dp_lcr(1:13),aic_d(1:13));
rFh = corr(dp_hcr(1:13),aic_d(1:13));

rJl = corr(dp_lcr(14:end),aic_d(14:end));
rJh = corr(dp_hcr(14:end),aic_d(14:end));

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
% => monkey F (N = 13)
sdFl = 1/sqrt(13 - 3); sdFh = 1/sqrt(13 - 3);
% => monkey J (N = 16)
sdJl = 1/sqrt(16 - 3); sdJh = 1/sqrt(16 - 3);

figure(20); set(gcf,'color','white');
hold on; hold all;
% ==> monkey F
errorbar([1,2], ...
         [zFli, zFhi], ...
         [zi(sdFl),zi(sdFh)] , ...
         'o','color','k');
% ==> monkey J
errorbar([3,4], ...
         [zJli, zJhi], ...
         [zi(sdJl),zi(sdJh)] , ...
         'ko');  

% figure(18); set(gcf,'color','white');
% hold on; hold all;
% % ==> monkey F
% errorbar([1,2], ...
%          [z(rFl), z(rFh)], ...
%          [sdFl,sdFh] , ...
%          'o','color','k');
% % ==> monkey J
% errorbar([3,4], ...
%          [z(rJl), z(rJh)], ...
%          [sdJl,sdJh] , ...
%          'ko');     
     
xlim([0,5]);
ylim([0,1.1]);