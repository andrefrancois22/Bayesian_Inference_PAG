clear all; close all; clc;
% ==> directories
dataPath     = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data';
drc = '../../data/';
% ==> directory containing PF curve fit functions
pfc_functions_dr = '/home/thomas/Desktop/Bayesian_Inference_PAG/simulation/F1/pfc_functions/';
addpath(pfc_functions_dr)

% ==> dynamic range analysis

% ==> dynamic range analysis types
dyn_types = {'popu','stim'};
% ==> estimate dynamic range per stimulus orientation? or per population?
dyn_type = dyn_types{1};

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

% ==> delta bias and delta uncertainty metrics
db = nan(29,2);
dp = nan(29,2);

% ==> PF curve fits
PFs_predPF_dynr_lo = cell(29,1);
PFs_predPF_dynr_hi = cell(29,1);

% ==> PF curve fits
PFs_predPF_dynr_lo_hcr = cell(29,1);
PFs_predPF_dynr_hi_hcr = cell(29,1);

figure(3); set(gcf,'color','white'); set(gcf,'Position',[86 282 921 668]);

% ==> what is the session
for iS = 1:29 

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

    % ==> choice
    cho = S.beh.choiceCat;
    % ==> context, contrast, orientation indicator variables
    ctx = S.exp.taskContext; ctr = S.exp.stimContrast; ori = S.exp.stimOriDeg;
    % x axis is orientation (the values differ for FN and JP!). Use unique values
    or = unique(ori)'; cx = unique(ctx)'; cr = unique(ctr)';
    
    % ==> split by dynamic range over entire population DVs
    if strcmp(dyn_type,'popu')
        % ==> dynamic range for session
        dynr = dynf(dvs); ldyn = (dynr(:,1) <  median(dynr)); hdyn = (dynr(:,1) >= median(dynr));
    end
    
    %for bh  (1 or 2) -> (-1 or 1)
        %for xc (1 or 2) -> (-1 or 1)
            %for rc (1 or 2) -> (min or max)
                %for sp (1 or 2) -> ldyn/hdyn
                
                    % idxs{bh}{xc}{rc}{sp} = (ctx == cx(xc)) & (ctr == cr)
                
                %end
            %end
        %end
    %end
    
    % ==> cw context (congruent choices)
    % ==> cw context & lo contrast
    cw_lcr_cng = (ctx == 1) & (ctr == min(ctr))   & (cho == 1); 
    % ==> ccw context & lo contrast
    ccw_lcr_cng = (ctx == -1) & (ctr == min(ctr)) & (cho == -1);      
    % ==> cw context (incongruent choices)
    % ==> cw context & lo contrast
    cw_lcr_icg = (ctx == 1) & (ctr == min(ctr))   & (cho == -1);           
    % ==> ccw context (incongruent choices)
    % ==> ccw context & lo contrast
    ccw_lcr_icg = (ctx == -1) & (ctr == min(ctr)) & (cho == 1);   
        
    % ==> cw context & hi contrast
    cw_hcr_cng = (ctx == 1) & (ctr == max(ctr))   & (cho == 1);     
    % ==> ccw context (congruent choices)
    % ==> ccw context & hi contrast
    ccw_hcr_cng = (ctx == -1) & (ctr == max(ctr)) & (cho == -1);      
    % ==> cw context (incongruent choices)
    % ==> cw context & hi contrast
    cw_hcr_icg = (ctx == 1) & (ctr == max(ctr))   & (cho == -1);              
    % ==> ccw context (incongruent choices)
    % ==> ccw context & hi contrast
    ccw_hcr_icg = (ctx == -1) & (ctr == max(ctr)) & (cho == 1);       
                       
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
        % ==> split by dynamic range conditioned on identical stimulus values
        if strcmp(dyn_type,'stim')
            % ==> dynamic range for session % => indexing variables - dynamic range
            dv = dvs((ctx == 1) & (ctr == min(ctr)) & ori==or(th),:);
            dynr = dynf(dv); ldyn=(dynf(dvs) <  median(dynr)); hdyn=(dynf(dvs) >= median(dynr));                
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
            dynr = dynf(dv); ldyn=(dynf(dvs) <  median(dynr)); hdyn=(dynf(dvs) >= median(dynr));             
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
            dynr = dynf(dv); ldyn=(dynf(dvs) <  median(dynr)); hdyn=(dynf(dvs) >= median(dynr));                         
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
            dynr = dynf(dv); ldyn=(dynf(dvs) <  median(dynr)); hdyn=(dynf(dvs) >= median(dynr));                       
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
    tn_ccw_lo_p = tn_ccw_lo ./ sum(tn_ccw_lo + tn_cw_lo); 
    tn_cw_lo_p =  tn_cw_lo  ./ sum(tn_ccw_lo + tn_cw_lo); 
    % ==> high dynamic range
    tn_ccw_hi_p = tn_ccw_hi ./ sum(tn_ccw_hi + tn_cw_hi);
    tn_cw_hi_p =  tn_cw_hi  ./ sum(tn_ccw_hi + tn_cw_hi);  
    
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
    
    % ==> delta bias (this is a different computation from original)
    db(iS,1) = (p_ldyn_lcr(5) - p_ldyn_lcr(4)) - (p_hdyn_lcr(5) - p_hdyn_lcr(4));
    %paramEst_M1_dynr_lo(5) - paramEst_M1_dynr_hi(5);    
    % ==> delta perceptual uncertainty
    dp(iS,1)   = p_ldyn_lcr(3) - p_hdyn_lcr(3);
    
    % ==> delta bias (this is a different computation from original)
    db(iS,2) = (p_ldyn_hcr(5) - p_ldyn_hcr(4)) - (p_hdyn_hcr(5) - p_hdyn_hcr(4));
    %paramEst_M1_dynr_lo_hcr(5) - paramEst_M1_dynr_hi)_hcr(5);    
    % ==> delta perceptual uncertainty
    dp(iS,2)   = p_ldyn_hcr(3) - p_hdyn_hcr(3);    
    
    % ==> save PF curve fits
    PFs_predPF_dynr_lo{iS} = PF_ldyn_lcr;
    PFs_predPF_dynr_hi{iS} = PF_hdyn_lcr;  
    % ==> save PF curve fits
    PFs_predPF_dynr_lo_hcr{iS} = PF_ldyn_hcr;
    PFs_predPF_dynr_hi_hcr{iS} = PF_hdyn_hcr;      
    
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


% ==> Average delta bias by animal and by contrast 
% ==> correlation plots (delta bias and delta AIC) 
[~,~,~,~,~,~,~] = F4EFG(db, aic_d,'bias');

% ==> average delta perceptual uncertainty
% ==> correlation plots (delta perceptual uncertainty and delta AIC
[~,~,~,~,~,~,~] = F4EFG(dp, aic_d,'uncertainty');

% ==> correlations between db and dp
F5(db, dp);

% ==> example plot
figure(); set(gcf,'color','white');
scatter(db(14:29,1), dp(14:29,1),75,'ko','filled');
xlabel('Difference in decision bias');
ylabel('Difference in slope');
axis square;

% % ==> ANCOVA
% dbz  = [db_lcr_z, db_hcr_z];
% dpz  = [dp_lcr_z, dp_hcr_z];
% aicz = [aic_d_z, aic_d_z];
% % => histogram to get grouping variable.
% [~,~,aicgv] = histcounts(aicz,8);
% % ==> run ANCOVA
% aoctool(dbz,dpz,aicgv)

