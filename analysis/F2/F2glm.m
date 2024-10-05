clear all; close all; clc;
% ==> directories
dataPath     = strcat('../../data/pfc_data/');
drc = '../../data/';

% ==> all sessions
Sn = 1:29;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ==> load dvCatperf data (for comparison to Delta AIC)
dvCatPerf = load([drc,'dvCatPerf.mat'],'dvCatPerf');
dvCatPerf = dvCatPerf.dvCatPerf;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~~~~~ comparing m0 with only orientation to m1 with orientation and signed peak
% => model 1 index for IVs
oidx_m1 = [1,3];
% => model 0 index for IVs
oidx_m0 = 1; 
% ==> AIC K params
% => m1 k = (2 x beta + beta0 (intercept)) x 4 models = 3 x 4 = 12
k1 = 12; 
% => m0 k = 4 x 2 = 8
k0 = 8;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> for storing 4 model design matrices
X = cell([29,2,2]);
% ==> y (categorical choice outcome)
Y = cell([29,2,2]);

% ==> what is the session
for iS = Sn

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
    or = sort(unique(ori),'ascend')'; cx = sort(unique(ctx),'ascend')'; cr = sort(unique(ctr),'ascend')';
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~ Cat DV (signed) peak ~~~~~~~~~~~~~~~~~~~~~~~ %
    % ==> get max (maximum is either the positive peak if it's highest, 
    % or abs of negative peak if it's lowest)
    mx = [max(dvs, [], 2), abs(min(dvs, [], 2))];
    mxv = max(mx, [], 2);
    % => signed
    mxv(mx(:,1) <= mx(:,2)) = -mxv(mx(:,1) <= mx(:,2));
    % ==> alpha is DV peak
    alpha = [sign(mxv), mxv];    
    
    for xc = 1:length(cx)
        for rc = 1:length(cr)
            % ==> index context xc, and contrast rc
            idx = (ctx == cx(xc)) & (ctr == cr(rc));
            % ==> store IVs
            X{iS,xc,rc} = [ori(idx), alpha(idx,:)];
            % ==> choices (turn '-1' to '0')
            y = cho(idx); y(y == -1) = 0;
            % ==> behavioral choice outcome for glm
            Y{iS,xc,rc} = y;
        end
    end
    % ==> update
    fprintf('Stored model data for session %d of %d...\n',iS,29)
end

%% 
clc;

% ==> log-likelihood ratio
llr = nan(1,29);

% ==> AIC delta (AIC_0 - AIC_1), where AIC_1 is model 1 (our preferred model)
aic_d = nan(29,1);
% ==> probability of the ith model
aic_p = nan(29,1);
% ==> BIC delta
bic_d = nan(29,1);

aic_m1 = cell([29,1]);
aic_m0 = cell([29,1]);

% ==> betas (model 1)
b1 = cell([29,2,2]);
% ==> y^hat (model 1)
yh1 = cell([29,2,2]);

% ==> betas (model 0)
b0 = cell([29,2,2]);
% ==> y^hat (model 0)
yh0 = cell([29,2,2]);

% ==> log likelihoods
m1ll = cell([29,2,2]);
m0ll = cell([29,2,2]);

% ==> models
m1 = cell([29,1]);
m0 = cell([29,1]);

% ==> log-likelihoods
ll1 = cell([29,1]);
ll0 = cell([29,1]);

for iS = Sn
    for xc = 1:length(cx)
        for rc = 1:length(cr)

            % ==> logistic regression (ori and alpha)
            [b1{iS,xc,rc},~,~] = glmfit(X{iS,xc,rc}(:,oidx_m1), Y{iS,xc,rc},'binomial','link','logit');    
            % => yhat (model 1)
            yh1{iS,xc,rc} = glmval(b1{iS,xc,rc}, X{iS,xc,rc}(:,oidx_m1),'logit');
            
            % ==> Null model (just orientation)
            [b0{iS,xc,rc},~,~] = glmfit(X{iS,xc,rc}(:,oidx_m0), Y{iS,xc,rc},'binomial','link','logit');    
            % => yhat (model 0)
            yh0{iS,xc,rc} = glmval(b0{iS,xc,rc}, X{iS,xc,rc}(:,oidx_m0),'logit');
            
            % => for computing log likelihood (model 1)
            m1ll{iS,xc,rc} = log(yh1{iS,xc,rc}).*Y{iS,xc,rc} + log(1 - yh1{iS,xc,rc}).*(1 - Y{iS,xc,rc});
            % => for computing log likelihood (model 0)
            m0ll{iS,xc,rc} = log(yh0{iS,xc,rc}).*Y{iS,xc,rc} + log(1 - yh0{iS,xc,rc}).*(1 - Y{iS,xc,rc});

            % ==> aggregate for computing log likelihood (with a sum over all observations)
            m1{iS} = [m1{iS}; m1ll{iS,xc,rc}]; 
            m0{iS} = [m0{iS}; m0ll{iS,xc,rc}];  
        end   
    end
    
    % ==> log likelihoods (lls)
    ll1{iS} = sum(m1{iS});
    ll0{iS} = sum(m0{iS});
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ==> log (likelihood ratio l1/l0) (difference in log-likelihoods e.g. ll1 - ll0)
    llr(iS) = ll1{iS} - ll0{iS};
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % ==> AIC (k values are numbers of parameters for each model)
    aic_m1{iS} = -(2/length(m1{iS}))*ll1{iS} + (2*k1/length(m1{iS}));    
    aic_m0{iS} = -(2/length(m0{iS}))*ll0{iS} + (2*k0/length(m0{iS}));

    % ==> choose model with lowest AIC. Quantify this with a difference
    % => values above 0 indicate model 1 is better than null model
    aic_d(iS) = (aic_m0{iS} - aic_m1{iS}); 

end

% ==> correlation between dvCatPerf and log likelihood ratios
[rho_llr,  p_llr] = corr(dvCatPerf',llr','type','spearman');
% ==> correlation between dvCatPerf and AIC delta
[rho_aicd, p_aicd] = corr(dvCatPerf',aic_d,'type','spearman');
% ==> print results (sanity check - compare predictions to LDA choice predictivity)
fprintf('r between dvCatPerf and log-likelihood ratios r = %d (p = %d)\n', rho_llr, p_llr)
fprintf('r between dvCatPerf and AIC Delta r = %d (p = %d)\n', rho_aicd, p_aicd)

% ==> save AIC delta and llr
save([drc,'aic_d.mat'],'aic_d');
save([drc,'llr.mat'],'llr');

% => monkey colors
clrs = {[1,0.75,0],[0.15,0.75,0.5]};
% => monkey session ranges
rgs = {1:13,14:29};
% => monkeys
M = {'F','J'};
% ==> draw figure
fg = figure();  set(fg,'Color','w'); fg.Position = [280 363 1084 482];
% ==> plot LLRs
for m = 1:length(M)
    ax = subplot(1,2,1);
    % ==> plot log-likelihood ratio results
    scatter(llr(rgs{m}),  rgs{m}, 150,  'o', 'filled',  'markeredgecolor', 'w', 'markerfacecolor', clrs{m});
    xlabel('log-likelihood ratio'); ylabel('Session'); xlim([-20,max(llr)]); ax.YTickLabel = [];
    hold on; hold all; plot(zeros(1,30),0:29,'k--');
    % ==> plot AIC Delta
    ax2 = subplot(1,2,2);
    scatter(aic_d(rgs{m}),  rgs{m},  150, 'o', 'filled', 'markeredgecolor', 'w', 'markerfacecolor', clrs{m});   
    xlabel('AIC \Delta'); ylabel('Session'); ax.YTickLabel = [];
    hold on; hold all; plot(zeros(1,30),0:29,'k--');    
end