clear all; close all; clc;
% ==> directories
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '../data/';

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ==> load dvCatperf data
dvCatPerf = load([drc,'dvCatPerf.mat'],'dvCatPerf');
dvCatPerf = dvCatPerf.dvCatPerf;
% => sort and return index
[~,idx] = sort(dvCatPerf);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> for storing 4 model design matrices
% => context 1 & contrast 1
cx1_cr1_X = cell([29,1]);
% => context 1 & contrast 2
cx1_cr2_X = cell([29,1]);

% => context 2 & contrast 1
cx2_cr1_X = cell([29,1]);
% => context 2 & contrast 2
cx2_cr2_X = cell([29,1]);

% ==> y (categorical choice outcome)
Y11 = cell([29,1]);
Y12 = cell([29,1]);
Y21 = cell([29,1]);
Y22 = cell([29,1]);

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
    % ==> behavioral choice outcome is y
    y11 = cho(idx11);
    y11(y11 == -1) = 0;    
    % => store
    Y11{iSfln} = y11;
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % => index context 1, and contrast 2 (hi contrast)
    idx12 = (ctx == 1) & (ctr == cr(2));
    % 2) ==> dummy coding orientation        
       
    % => context 1 & contrast 1
    cx1_cr2_X{iSfln} = [ori(idx12), dynran(idx12)];    
    % ==> behavioral choice outcome is y
    y12 = cho(idx12);
    y12(y12 == -1) = 0;    
    % => store
    Y12{iSfln} = y12;
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % => index context 2, and contrast 1 (low contrast)
    idx21 = (ctx == -1) & (ctr == cr(1));
    % 2) ==> dummy coding orientation        
       
    % => context 1 & contrast 1
    cx2_cr1_X{iSfln} = [ori(idx21), dynran(idx21)];
    % ==> behavioral choice outcome is y
    y21 = cho(idx21);
    y21(y21 == -1) = 0;   
    % => store
    Y21{iSfln} = y21;
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
    % => store
    Y22{iSfln} = y22;
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
    fprintf('Stored model data for session %d of %d...\n',iSfln,29)
end

%% 
close all; clc;

% => orientation independent variables index
oidx = 1;

% ==> log-likelihood ratio
llr = nan(1,29);

% ==> AIC delta (AIC_0 - AIC_1), where AIC_1 is model 1 (our preferred model)
aic_d = nan(1,29);
% ==> probability of the ith model
aic_p = nan(1,29);
% ==> BIC delta
bic_d = nan(1,29);

for iS = 1:29

    % ====> condition 11 (context 1 - low contrast)
    % ==> logistic regression (ori and dynran)
    [b_m1_11,~,~] = glmfit(cx1_cr1_X{iS}, Y11{iS},'binomial','link','logit');    
    % => yhat
    yh_m1_11 = glmval(b_m1_11, cx1_cr1_X{iS},'logit');
    % ==> Null model (just orientation)
    [b_m0_11,~,~] = glmfit(cx1_cr1_X{iS}(:,oidx), Y11{iS},'binomial','link','logit');    
    % => yhat
    yh_m0_11 = glmval(b_m0_11, cx1_cr1_X{iS}(:,oidx),'logit');
    % => for computing log likelihood (model 1)
    m1ll_11 = log(yh_m1_11).*Y11{iS} + log(1 - yh_m1_11).*(1 - Y11{iS});
    % => for computing log likelihood (model 0)
    m0ll_11 = log(yh_m0_11).*Y11{iS} + log(1 - yh_m0_11).*(1 - Y11{iS});


    % ====> condition 12 (context 1 - hi contrast)
    % ==> logistic regression (ori and dynran)
    [b_m1_12,~,~] = glmfit(cx1_cr2_X{iS}, Y12{iS},'binomial','link','logit');    
    % => yhat
    yh_m1_12 = glmval(b_m1_12, cx1_cr2_X{iS},'logit');
    % ==> Null model (just orientation)
    [b_m0_12,~,~] = glmfit(cx1_cr2_X{iS}(:,oidx), Y12{iS},'binomial','link','logit');    
    % => yhat
    yh_m0_12 = glmval(b_m0_12, cx1_cr2_X{iS}(:,oidx),'logit');
    % => for computing log likelihood (model 1)
    m1ll_12 = log(yh_m1_12).*Y12{iS} + log(1 - yh_m1_12).*(1 - Y12{iS});
    % => for computing log likelihood (model 0)
    m0ll_12 = log(yh_m0_12).*Y12{iS} + log(1 - yh_m0_12).*(1 - Y12{iS});


    % ====> condition 21 (context 2 - low contrast)
    % ==> logistic regression (ori and dynran)
    [b_m1_21,~,~] = glmfit(cx2_cr1_X{iS}, Y21{iS},'binomial','link','logit');    
    % => yhat
    yh_m1_21 = glmval(b_m1_21, cx2_cr1_X{iS},'logit');
    % ==> Null model (just orientation)
    [b_m0_21,~,~] = glmfit(cx2_cr1_X{iS}(:,oidx), Y21{iS},'binomial','link','logit');    
    % => yhat
    yh_m0_21 = glmval(b_m0_21, cx2_cr1_X{iS}(:,oidx),'logit');
    % => for computing log likelihood (model 1)
    m1ll_21 = log(yh_m1_21).*Y21{iS} + log(1 - yh_m1_21).*(1 - Y21{iS});
    % => for computing log likelihood (model 0)
    m0ll_21 = log(yh_m0_21).*Y21{iS} + log(1 - yh_m0_21).*(1 - Y21{iS});

    
    % ====> condition 22 (context 2 - hi contrast)
    % ==> logistic regression (ori and dynran)
    [b_m1_22,~,~] = glmfit(cx2_cr2_X{iS}, Y22{iS},'binomial','link','logit');    
    % => yhat
    yh_m1_22 = glmval(b_m1_22, cx2_cr2_X{iS},'logit');
    % ==> Null model (just orientation)
    [b_m0_22,~,~] = glmfit(cx2_cr2_X{iS}(:,oidx), Y22{iS},'binomial','link','logit');    
    % => yhat
    yh_m0_22 = glmval(b_m0_22, cx2_cr2_X{iS}(:,oidx),'logit');
    % => for computing log likelihood (model 1)
    m1ll_22 = log(yh_m1_22).*Y22{iS} + log(1 - yh_m1_22).*(1 - Y22{iS});
    % => for computing log likelihood (model 0)
    m0ll_22 = log(yh_m0_22).*Y22{iS} + log(1 - yh_m0_22).*(1 - Y22{iS});    
    
    % ==> aggregate for computing log likelihood (with a sum over all observations)
    m1 = [m1ll_11; m1ll_21; m1ll_12; m1ll_22];
    m0 = [m0ll_11; m0ll_21; m0ll_12; m0ll_22];
    % ==> log likelihoods (lls)
    ll1 = sum(m1);
    ll0 = sum(m0);
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ==> log (likelihood ratio l1/l0) (difference in log likelihoods e.g. ll1 - ll0)
    llr(iS) = ll1 - ll0;
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % ==> AIC
    % => m1 (k = 2 betas + beta0 (intercept) + epsilon noise) == 4
    k1 = 4;
    % => m0 (k = 1 beta + beta0 (intercept) + epsilon noise) == 3
    k0 = 3;
    aic_m1 = -2/length(m1)*ll1 + 2*k1/length(m1);    
    aic_m0 = -2/length(m0)*ll0 + 2*k0/length(m0);
    
    % ==> choose model with lowest AIC. Quantify this with a difference
    % => values above 0 indicate model 1 is better than null model
    aic_d(iS) = (aic_m0 - aic_m1); 
    %aic_p(iS) = exp( (aic_m0 - aic_m1) * 0.5 );
    
    % ==> BIC
    bic_m1 = -2*ll1 + log(length(m1))*k1;
    bic_m0 = -2*ll0 + log(length(m0))*k0;    
    % ==> choose model with lowest BIC. Quantify this with a difference
    % => values above 0 indicate model 1 is better than null model
    bic_d(iS) = (bic_m0 - bic_m1); % ==> want bic_m1 to be less than bic_m0
    
    
end

% ==> correlation between dvCatPerf and log likelihood ratios
[rho_llr,  p_llr] = corr(dvCatPerf(idx)',llr(idx)','type','spearman'); 
% ==> correlation between dvCatPerf and AIC delta
[rho_aicd, p_aicd] = corr(dvCatPerf(idx)',aic_d(idx)','type','spearman');
% ==> correlation between dvCatPerf and BIC delta
[rho_bicd, p_bicd] = corr(dvCatPerf(idx)',bic_d(idx)','type','spearman');

fprintf('r between dvCatPerf and log-likelihood ratios r = %d (p = %d)\n', rho_llr, p_llr)
fprintf('r between dvCatPerf and AIC Delta r = %d (p = %d)\n', rho_aicd, p_aicd)
fprintf('r between dvCatPerf and BIC Delta r = %d (p = %d)\n', rho_bicd, p_bicd)

fg = figure();  set(fg,'Color','w'); fg.Position = [675 241 512 721];
ax1 = subplot(3,2,1);
scatter(llr, 1:29, 'filled'); 
xlabel('log-likelihood ratio'); ylabel('Session'); xlim([-20,max(llr)]); ax1.YTickLabel = [];
hold on; hold all; plot(zeros(1,30),0:29,'k--');
subplot(3,2,2);
scatter(dvCatPerf(idx)',llr(idx)', 'filled');
xlabel('Choice Predictivity'); ylabel('log-likelihood ratio');
ax3 = subplot(3,2,3);
scatter(aic_d, 1:29, 'filled');
xlabel('AIC \Delta'); ylabel('Session');  xlim([-0.02,max(aic_d)]); ax3.YTickLabel = [];
hold on; hold all; plot(zeros(1,30),0:29,'k--');
subplot(3,2,4);
scatter(dvCatPerf(idx)',aic_d(idx)', 'filled'); 
xlabel('Choice Predictivity'); ylabel('AIC \Delta');
ax5 = subplot(3,2,5);
scatter(bic_d, 1:29, 'filled'); 
xlabel('BIC \Delta'); ylabel('Session');  xlim([-20,max(bic_d)]); ax5.YTickLabel = [];
hold on; hold all; plot(zeros(1,30),0:29,'k--');
subplot(3,2,6);
scatter(dvCatPerf(idx)',bic_d(idx)', 'filled'); 
xlabel('Choice Predictivity'); ylabel('BIC \Delta');


fg = figure();  set(fg,'Color','w'); fg.Position = [670 219 717 731];
ax1 = subplot(2,2,1);
scatter(llr, 1:29, 'filled'); 
xlabel('log-likelihood ratio'); ylabel('Session'); xlim([-20,max(llr)]); ax1.YTickLabel = [];
hold on; hold all; plot(zeros(1,30),0:29,'k--');
subplot(2,2,2);
scatter(dvCatPerf(idx)',llr(idx)', 'filled');
xlabel('Choice Predictivity'); ylabel('log-likelihood ratio'); title(['r = ',num2str(rho_llr),', p = ', num2str(p_llr)]);
ax3 = subplot(2,2,3);
scatter(aic_d, 1:29, 'filled');
xlabel('AIC \Delta'); ylabel('Session');  xlim([-0.02,max(aic_d)]); ax3.YTickLabel = [];
hold on; hold all; plot(zeros(1,30),0:29,'k--');
subplot(2,2,4);
scatter(dvCatPerf(idx)',aic_d(idx)', 'filled');  title(['r = ',num2str(rho_aicd),', p = ', num2str(p_aicd)]); 
xlabel('Choice Predictivity'); ylabel('AIC \Delta');

