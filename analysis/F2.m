clear all; close all; clc;
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/goris_data/data/';

% => pairwise combinations 
mrx = nan(29,2,2,7);

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

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % ~~~~~~~~~~~~~~~~~~~~~~~~ Cat DV (signed) peak ~~~~~~~~~~~~~~~~~~~~~~~ %
    % ==> get max (maximum is either the positive peak if it's highest, or abs of negative peak if it's lowest)
    mx = [max(dvs, [], 2), abs(min(dvs, [], 2))];
    mxv = max(mx, [], 2);
    % => signed
    mxv(mx(:,1) <= mx(:,2)) = -mxv(mx(:,1) <= mx(:,2));
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    % x axis is orientation (the values differ for FN and JP!). Use unique values
    or = unique(ori)'; cx = unique(ctx)'; cr = unique(ctr)';
    % => store averages for each stimulus condition
    % => loop through two context conditions
    for j = 1:length(cx)
        % => loop through two contrast conditions
        for k = 1:length(cr)
            % => loop through orientations
            for i = 1:length(or)
                mrx(iSfln,j,k,i) = mean( mxv( ori==or(i) & ctr==cr(k) & ctx==cx(j) ) );
            end                      
        end                    
    end      
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % ~~~~~~~~~~~~~~~~~~~~~~~~ Plot results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    clrs = {'b','r'}; mrks = {'o:','*-'};
    figure(1); set(gcf,'color','white'); 
    subplot(5,6,iSfln); axis square;
    for j = 1:length(cx)
        for k = 1:length(cr)
            plot(or,squeeze(mrx(iSfln,j,k,:))',[clrs{j},mrks{k}]); hold on; hold all;        
        end
    end 
    title(['Session: ',num2str(iSfln)]);
    ylabel('Peak Average (signed)')
    xlabel('Orientation')
    drawnow;
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    
    % => console updates...
    fprintf('Computed results for session %d of %d... \n',iSfln,29)
end

%% Fit bias and slope to DV peaks

close all; clc;
% ==> bias and perceptual uncertainty vectors
bv  = nan(29,2);
puv = nan(29,2); 

% ==> for storing just slopes
slope = nan(29,2);
% => track loss for sanity check
lss = nan(29,2);

for iSfln =  1:29 

clrs = {'b','r'}; mrks = {'o:','*-'};
figure(2); set(gcf,'color','white'); 
clf; axis square;

% => k is contrast
for k = 1:2    
    
    % ==> fit linear model with identical slopes and different intercepts
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ==> ys and x ir orientations
    y = squeeze(mrx(iSfln,:,k,:));
    x = or;
    % ==> fit linear models with identical slopes and different intercepts to
    % => initial parameter settings
    p0 = [0,0,0];
    % ==> two lines - 
    yh1 = @(p) (p(1)*x + p(2));
    yh2 = @(p) (p(1)*x + p(2) + p(3));
    % ==> MSE objective for both lines
    objective = @(p) sum( (yh1(p) - y(1,:)).^2 + (yh2(p) - y(2,:)).^2 );
    % ==> overkill - could probably use normal equations optim params...
    p_opt = fmincon(objective,p0);
    % => inspect loss
    lss(iSfln,k) = objective(p_opt);        

    % ==> the bias is the vertical offset (optimal parameter 3 p(3))
    bv(iSfln,k)  = p_opt(3);
    % ==> uncertainty is slope
    puv(iSfln,k) = 1 / p_opt(1); %**** of course it's the inverse! ==> shallower slope is ~larger~ likelihood width      
    % ==> store slopes
    slope(iSfln,k) = p_opt(1);
    
    % ==> console update
    fprintf('completed linear model fit for context %d and contrast %k...\n',j,k)           
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
    % ==> for visibility to inspect curve fits
    subplot(1,2,k)
    for j = 1:length(cx)  
        plot(or,squeeze(mrx(iSfln,j,k,:))',[clrs{j},mrks{k}],'linewidth',1); hold on; hold all;   
        if j==1
            % ==> plot line fit          
            plot(or,yh1(p_opt),['k','h:']);
        elseif j==2
            % ==> plot line fit  
            plot(or,yh2(p_opt),['k','>:']);   
        end
    end
    title('Peak Average (signed) by stimulus condition and orientation')    
    title(['Session: ',num2str(iSfln)]);
    ylabel('Peak Average (signed)');
    xlabel('Orientation');    
    ylim([-2,2]);    
end
%       pause(1)
end

% ==> display loss
fprintf('total fit loss MSE = %d\n...',mean(lss(:)))
b = bv(~isnan(bv(:)));
p = puv(~isnan(puv(:)));

% F: 1-13; 
bF = bv(1:13,:);  bF = bF(:);     
pF = puv(1:13,:); pF = pF(:);
% ==> r
[rF, pvF] = corr(real(log10(pF)),bF,'type','pearson');
% JP: 14-29
bJ = bv(14:end,:);  bJ = bJ(:);     
pJ = puv(14:end,:); pJ = pJ(:);
% ==> r
[rJ,pvJ] = corr(real(log10(pJ)),bJ,'type','pearson');
clc;
fprintf('Rank correlation for Friederich is r = %d, and for Jean-Paul is r = %d...\n',rF,rJ)


% F: 1-13; JP: 14-29
% ==> plot and rank correlation
figure(); set(gcf,'color','white');
% ==> spearman correlation
semilogx(puv(1:13,1),bv(1:13,1),'o', 'markeredgecolor', [0.15,0.75,0.5], 'markerfacecolor', [1,1,1], 'linewidth', 2, 'markersize',10); axis square;
hold on; hold all;
semilogx(puv(1:13,2),bv(1:13,2),'o', 'markeredgecolor', [0.15,0.75,0.5], 'markerfacecolor', [0.15,0.75,0.5], 'linewidth', 2, 'markersize',10); axis square; 
hold on; hold all;
semilogx(puv(14:end,1),bv(14:end,1),'o', 'markeredgecolor', [1,0.5,0], 'markerfacecolor', [1,1,1], 'linewidth', 2, 'markersize',10); axis square; 
hold on; hold all;
semilogx(puv(14:end,2),bv(14:end,2),'o', 'markeredgecolor', [1,0.5,0], 'markerfacecolor', [1,0.5,0], 'linewidth', 2, 'markersize',10); axis square;
xlabel('log_{10} (perceptual uncertainty (1 / \beta_1))','fontsize',20)
ylabel('bias (\beta_0)','fontsize',20)
title(['r = ',num2str(corr(real(log10(p)),b,'type','spearman'))]);
legend('Monkey F low contrast', 'Monkey F high contrast', ...
       'Monkey J low contrast', 'Monkey J high contrast')

% F: 1-13; JP: 14-29
% ==> plot and rank correlation
fg = figure(); set(fg,'color','white'); fg.Position = [675 309 1221 653];
subplot(1,2,1); 
% ==> spearman correlation
semilogx(puv(1:13,1),bv(1:13,1),'o', 'markeredgecolor', [0.15,0.75,0.5], 'markerfacecolor', [1,1,1], 'linewidth', 2, 'markersize',10); axis square;
hold on; hold all;
semilogx(puv(1:13,2),bv(1:13,2),'o', 'markeredgecolor', [0.15,0.75,0.5], 'markerfacecolor', [0.15,0.75,0.5], 'linewidth', 2, 'markersize',10); axis square; 
xlabel('log_{10} (perceptual uncertainty (1 / \beta_1))','fontsize',20)
ylabel('bias (\beta_0)','fontsize',20)
title(['r = ',num2str(rF), ' p = ',num2str(pvF)]);
legend('low contrast','High contrast')

subplot(1,2,2); 
semilogx(puv(14:end,1),bv(14:end,1),'o', 'markeredgecolor', [1,0.5,0], 'markerfacecolor', [1,1,1], 'linewidth', 2, 'markersize',10); axis square; 
hold on; hold all;
semilogx(puv(14:end,2),bv(14:end,2),'o', 'markeredgecolor', [1,0.5,0], 'markerfacecolor', [1,0.5,0], 'linewidth', 2, 'markersize',10); axis square;
xlabel('log_{10} (perceptual uncertainty (1 / \beta_1))','fontsize',20)
ylabel('bias (\beta_0)','fontsize',20)
title(['r = ',num2str(rJ), ' p = ',num2str(pvJ)]);
legend('low contrast','High contrast')


% => rank correlation with log10 - (pearson)

% ==> correlations by contrast conditions
[rFlo, pFlo] = corr(real(log10(puv(1:13,1))),  bv(1:13,1),'type','pearson');
[rJlo, pJlo] = corr(real(log10(puv(14:end,1))),bv(14:end,1),'type','pearson');

[rFhi, pFhi] = corr(real(log10(puv(1:13,2))),  bv(1:13,2),'type','pearson');
[rJhi, pJhi] = corr(real(log10(puv(14:end,2))),bv(14:end,2),'type','pearson');

figure(); set(gcf,'color','white');
subplot(2,2,1);
semilogx(puv(14:end,1),bv(14:end,1),'o', 'markeredgecolor', [1,0.5,0], 'markerfacecolor', [1,1,1], 'linewidth', 2, 'markersize',10); axis square; 
hold on; hold all;
xlabel('log_{10} (perceptual uncertainty (1 / \beta_1))','fontsize',10)
ylabel('bias (\beta_0)','fontsize',10)
title(['r = ',num2str(rJlo), ' p = ',num2str(pJlo)]);
subplot(2,2,2);
semilogx(puv(14:end,2),bv(14:end,2),'o', 'markeredgecolor', [1,0.5,0], 'markerfacecolor', [1,0.5,0], 'linewidth', 2, 'markersize',10); axis square;
xlabel('log_{10} (perceptual uncertainty (1 / \beta_1))','fontsize',10)
ylabel('bias (\beta_0)','fontsize',10)
title(['r = ',num2str(rJhi), ' p = ',num2str(pJhi)]);

subplot(2,2,3);
semilogx(puv(1:13,1),bv(1:13,1),'o', 'markeredgecolor', [0.15,0.75,0.5], 'markerfacecolor', [1,1,1], 'linewidth', 2, 'markersize',10); axis square;
hold on; hold all;
xlabel('log_{10} (perceptual uncertainty (1 / \beta_1))','fontsize',10)
ylabel('bias (\beta_0)','fontsize',10)
title(['r = ',num2str(rFlo), ' p = ',num2str(pFlo)]);
subplot(2,2,4);
semilogx(puv(1:13,2),bv(1:13,2),'o', 'markeredgecolor', [0.15,0.75,0.5], 'markerfacecolor', [0.15,0.75,0.5], 'linewidth', 2, 'markersize',10); axis square; 
xlabel('log_{10} (perceptual uncertainty (1 / \beta_1))','fontsize',10)
ylabel('bias (\beta_0)','fontsize',10)
title(['r = ',num2str(rFhi), ' p = ',num2str(pFhi)]);

%% Statistics for slope and bias differences

clc;
% ==> high contrast slopes minus low contrast slopes
sdiff = slope(:,2) - slope(:,1);

% ==>  signrank Wilcoxon signed rank test for zero median for F
[pF,hF] = signrank(sdiff(1:13));
fprintf('F slope difference Wilcoxon test p = %d, H = %d \n',pF, hF)
% ==>  signrank Wilcoxon signed rank test for zero median for JP
[pJ,hJ] = signrank(sdiff(14:end));
fprintf('JP slope difference Wilcoxon test p = %d, H = %d \n',pJ, hJ)
fprintf('...\n')

% ==> median F
mdF = median(sdiff(1:13));
% ==> median JP
mdJ = median(sdiff(14:end));

% ==> bias difference
bdiff = bv(:,2) - bv(:,1);

% ==>  signrank Wilcoxon signed rank test for zero median for F
[pbF,hbF] = signrank(bdiff(1:13));
fprintf('F bias difference Wilcoxon test p = %d, H = %d \n',pbF,hbF)
% ==>  signrank Wilcoxon signed rank test for zero median for JP
[pbJ,hbJ] = signrank(bdiff(14:end));
fprintf('JP bias difference Wilcoxon test p = %d, H = %d \n',pbJ,hbJ)

% ==> median F
mdbF = median(bdiff(1:13));
% ==> median JP
mdbJ = median(bdiff(14:end));

%% ==> F2 - slope Hi v slope Lo, and bias Hi v bias Lo

close all; clc;

figure(); set(gcf,'color','white');
subplot(1,4,1); 
scatter(slope(14:end,1),slope(14:end,2), 28,'o','filled'); axis square; ylim([0,0.6]); xlim([0,0.6]);
hold on; hold all;
plot(linspace(-0.1,0.65,10),linspace(-0.1,0.65,10),'k--')
xlabel('Slope low contrast')
ylabel('Slope high contrast')
title('JP')
subplot(1,4,2); 
scatter(slope(1:13,1),slope(1:13,2), 28,'o','filled'); axis square; ylim([0,0.4]); xlim([0,0.4]);
hold on; hold all;
plot(linspace(-0.1,0.65,10),linspace(-0.1,0.65,10),'k--')
xlabel('Slope low contrast')
ylabel('Slope high contrast')
title('F')

subplot(1,4,3); 
scatter(bv(14:end,1),bv(14:end,2), 28,'o','filled'); axis square; ylim([-0.4,1.1]); xlim([-0.4,1.1]);
hold on; hold all;
plot(linspace(-0.4,1,10),linspace(-0.4,1,10),'k--')
xlabel('Bias low contrast')
ylabel('Bias high contrast')
title('JP')
subplot(1,4,4); 
scatter(bv(1:13,1),bv(1:13,2), 28,'o','filled'); axis square; ylim([-0.1,0.65]); xlim([-0.1,0.65]);
hold on; hold all;
plot(linspace(-0.1,0.65,10),linspace(-0.1,0.65,10),'k--')
xlabel('Bias low contrast')
ylabel('Bias high contrast')
title('F')


%% ==> plot by choice congruency

clear all; close all; clc;
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/goris_data/data/';

dvs_all = [];
cw_i = [];
cw_c = [];
ccw_i = [];
ccw_c = [];

% ==> for dvCatPerf analysis
cw_dyn_dfs = [];
ccw_dyn_dfs = [];

for iSfln = 1:29

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

% ==> using orientation
% idx_i_cw = ori < 0 & (cho == 1);
% idx_c_cw = ori > 0 & (cho == 1);
% 
% idx_c_ccw = ori < 0 & (cho == -1);
% idx_i_ccw = ori > 0 & (cho == -1);

% ==> using context
idx_i_cw =  (ctx == -1) & (cho == 1);
idx_c_cw =  (ctx ==  1) & (cho == 1);

idx_c_ccw = (ctx == -1) & (cho == -1);
idx_i_ccw = (ctx ==  1) & (cho == -1);

% % ==> using context and ambiguous orientation
% % => absolute values of orientation degrees
% oris = unique(abs(ori));
% 
% idx_i_cw =  (ctx == -1) & (cho == 1) & (abs(ori) <= oris(1));
% idx_c_cw =  (ctx ==  1) & (cho == 1) & (abs(ori) <= oris(1));
% 
% idx_c_ccw = (ctx == -1) & (cho == -1) & (abs(ori) <= oris(1));
% idx_i_ccw = (ctx ==  1) & (cho == -1) & (abs(ori) <= oris(1));

% => append
cw_i = [cw_i; idx_i_cw];    % ==> incongruent case (true ccw stimulus)
cw_c = [cw_c; idx_c_cw];
ccw_i = [ccw_i; idx_i_ccw]; % ==> incongruent case (true cw stimulus)
ccw_c = [ccw_c; idx_c_ccw];

% => append
dvs_all = [dvs_all; dvs];

% => update
fprintf('completed session %d of %d...\n', iSfln,29)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ==> check this...
% ==> could do the dvCat perf analysis here too.
cw_i_dnr = mean(dvs_all(logical(cw_i),:),1);     cw_i_dnr_d = max(cw_i_dnr) - cw_i_dnr(1);
cw_c_dnr = mean(dvs_all(logical(cw_c),:),1);     cw_c_dnr_d = max(cw_c_dnr) - cw_c_dnr(1);

ccw_i_dnr = mean(dvs_all(logical(ccw_i),:),1);  ccw_i_dnr_d = abs(min(ccw_i_dnr) - ccw_i_dnr(1));
ccw_c_dnr = mean(dvs_all(logical(ccw_c),:),1);  ccw_c_dnr_d = abs(min(ccw_c_dnr) - ccw_c_dnr(1));

% ==> dyn rng difference between congruent and incongruent (cw)
% =>  positive value indicates higher dynamic range for incongruent cases
cw_dyn_df = cw_i_dnr_d / ccw_c_dnr_d;

% ==> dyn rng difference between congruent and incongruent (ccw)
% =>  positive value indicates higher dynamic range for incongruent cases
ccw_dyn_df = ccw_i_dnr_d / cw_c_dnr_d;

% ==> store
cw_dyn_dfs  = [cw_dyn_dfs,   cw_dyn_df];
ccw_dyn_dfs = [ccw_dyn_dfs, ccw_dyn_df];
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end

timeSac = S.dec.timeSac;
% Set time boundaries
fitTimeSacBegin = -800;  % in ms, relative to saccade onset
fitTimeSacEnd   = -50;
[~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
[~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));

% time intervals
t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));

figure;  set(gcf,'color','white');
plot(t, mean(dvs_all(logical(cw_i),:),1),'b-','linewidth',2) %=> incongr
hold on; hold all;
plot(t, mean(dvs_all(logical(cw_c),:),1),'r.-','linewidth',5)
hold on; hold all;
plot(t, mean(dvs_all(logical(ccw_i),:),1),'r-','linewidth',2) %=> incongr
hold on; hold all;
plot(t, mean(dvs_all(logical(ccw_c),:),1),'b.-','linewidth',5)
xticks(timeSac);
xticklabels(num2cell(timeSac))
xlim([timeSac(sacBegin),timeSac(sacEnd)]);
% => horizontal line at zero
plot(timeSac,zeros(1,size(timeSac,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[-0.5,0.5],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[-0.5,0.5],'k:')
% ylim([-0.5,0.5]);
xlabel('Time');
ylabel('Signed Categorical DV');
title('Categorical Decision Variables (DVs)');
axis square;

% ==> dvCatPerf
dvCatPerf = load([drc,'dvCatPerf.mat'],'dvCatPerf');
dvCatPerf = dvCatPerf.dvCatPerf;
% ==> plot
figure(); set(gcf,'color','white');
subplot(1,2,1);
plot(dvCatPerf,cw_dyn_dfs,'ro','markerfacecolor','r','markeredgecolor','r'); axis square;
xlabel('dvCatPerf');
ylabel('dynamic range ratio (cw)')
title([ 'r = ', num2str(corr(dvCatPerf',cw_dyn_dfs','type','spearman'))]);
subplot(1,2,2);
plot(dvCatPerf,ccw_dyn_dfs,'bo','markerfacecolor','b','markeredgecolor','b'); axis square;
xlabel('dvCatPerf');
ylabel('dynamic range ratio (ccw)')
title([ 'r = ', num2str(corr(dvCatPerf',ccw_dyn_dfs','type','spearman'))]);

[r_cw,  p_cw]  = corr(dvCatPerf', cw_dyn_dfs', 'type','spearman')
[r_ccw, p_ccw] = corr(dvCatPerf', ccw_dyn_dfs','type','spearman')

figure(); set(gcf,'color','white');
% ==> plot all points across cw and ccw
plot([dvCatPerf';dvCatPerf'],[cw_dyn_dfs'; ccw_dyn_dfs'],'ko' ,'markerfacecolor','k','markeredgecolor','k'); axis square;
xlabel('dvCatPerf');
ylabel('dynamic range ratio (ccw)')
title([ 'r = ', num2str(corr([dvCatPerf';dvCatPerf'],[cw_dyn_dfs'; ccw_dyn_dfs'],'type','spearman'))]);

[r_all, p_all] = corr([dvCatPerf';dvCatPerf'],[cw_dyn_dfs'; ccw_dyn_dfs'],'type','spearman')

%% F3 ?

% => plot
figure(); set(gcf,'color','white');

% different alphas
alphas = {'start_height', 'dynamic_range', 'dynamic_range_multiply', ...
          'rise_speed', 'peak_height', 'rise_time'};
      
for j = 1:length(alphas)

    alpha_param = alphas{j};
    
    % => session
    for iSfln = 1:29

        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % => single trial alpha values
        av = load([drc,'trial_DV_traj_alphas_iS_',num2str(iSfln),'_',alpha_param,'.mat'],  'av'); 
        av = av.av;
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % ==> load session data (would be faster to just do this once for all sessions)
        if iSfln > 9
            load([dataPath,'/dataSet_',num2str(iSfln),'.mat']);
        elseif iSfln <= 9
            load([dataPath,'/dataSet_0',num2str(iSfln),'.mat']);
        end
        % ==> context
        ctx = S.exp.taskContext; 
        fprintf('loaded dataSet iS=%d...\n',iSfln)
        % => total number of trials
        N = size(ctx,1);

        % => signed alphas need to be unsigned for mean difference
        av = abs(av);

        % mean and SE
        y = mean(av(ctx==1)) - mean(av(ctx==-1));
        % => standard error (SE) with pooled standard deviation
        SE = ( sqrt( (var(av(ctx==1)) + var(av(ctx==-1)))/2) ) / sqrt(N);

        % ==> plot
        subplot(2,3,j); grid on;
        title(['alpha = ',alpha_param], 'Interpreter','none');
        axis square; hold on; hold all;
        ylabel('Session'); 
        yticks([0:30]); ylim([0,30]);
        errorbar(y,iSfln,SE,'r.-','horizontal')
        drawnow;

    end
end

%% F2B updated

clear all; close all; clc;
dataPath = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/goris_data/data/';

% => one way is to store the indexed DVs across sessions
dvs_cw_or1_r = [];
dvs_cw_or2_r = [];
dvs_cw_or3_r = [];
dvs_cw_or0_r = [];
dvs_ccw_or1_r = [];
dvs_ccw_or2_r = [];
dvs_ccw_or3_r = [];

% => DV model fits
dvs_cw_or1_m = [];
dvs_cw_or2_m = [];
dvs_cw_or3_m = [];
dvs_cw_or0_m = [];
dvs_ccw_or1_m = [];
dvs_ccw_or2_m = [];
dvs_ccw_or3_m = [];

rs = nan(29,7);

figure();
% ==> what is the session
for iSfln = 14:29%     %1:13 % 
    
%     alpha_param = 'dynamic_range';

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
    
%     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     % => single trial alpha values
%     av = load([drc,'trial_DV_traj_alphas_iS_',num2str(iSfln),'_',alpha_param,'.mat'],  'av'); 
%     av = av.av;
%     % => indices (4 columns: 1=>idx_ccw_lo, 2=>idx_ccw_hi, 3=>idx_cw_lo, 4=> idx_cw_hi)
%     idxm = load([drc,'trial_medSplit_idxVec_',num2str(iSfln),'_',alpha_param,'.mat'],  'idxm'); 
%     idxm = idxm.idxm;
%     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    ctx = S.exp.taskContext; 
    ctr = S.exp.stimContrast; 
    ori = S.exp.stimOriDeg;

    % => contrast values (changes depending on JP of FN!)
    ctrv = unique(ctr);
    % => unique orientations 
    oriv = unique(ori); assert(length(oriv) == 7);

    % =-> actual Cat DVs
    dvs_real = S.dec.dvCatAll;

    % ==> subplots
    subplot(5,6,iSfln); hold on; hold all;

    % ==> high contrast, all orientations
    % => 1.1 degrees (JP)
    cw_or1 = (ori == oriv(end-2))           & (ctr==ctrv(2));
    % => 2.2 degrees (JP)
    cw_or2 = (ori == oriv(end-1))           & (ctr==ctrv(2));
    % => 3.3 degrees (JP)
    cw_or3 = (ori == oriv(end))             & (ctr==ctrv(2));
    % => 0 degrees
    cw_or0 = (ori == oriv(end-3))           & (ctr==ctrv(2));
    
    % => 11.1 degrees (JP)    
    ccw_or1 = (ori == oriv(3))              & (ctr==ctrv(2));
    % => -2.2 degrees (JP)
    ccw_or2 = (ori == oriv(2))              & (ctr==ctrv(2));
    % => -3.3 degrees (JP)
    ccw_or3 = (ori == oriv(1))              & (ctr==ctrv(2));
            
    % ==> plot 
    plot(1:size(dvs_real,2),mean(dvs_real(cw_or1,:),1),'r-'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(cw_or2,:),1),'r-'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(cw_or3,:),1),'r:'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(cw_or0,:),1),'r:'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(ccw_or1,:),1),'b-'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(ccw_or2,:),1),'b-'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(ccw_or3,:),1),'b:'); hold on; hold all;
    
    ylim([-1,1]);
    drawnow;
    
    % ==> for F1 plot 
    % => cw
    dvs_cw_or1_r = [dvs_cw_or1_r; dvs_real(cw_or1,:)];    
    dvs_cw_or1_m = [dvs_cw_or1_m;      dvs(cw_or1,:)];
    
    dvs_cw_or2_r = [dvs_cw_or2_r; dvs_real(cw_or2,:)];    
    dvs_cw_or2_m = [dvs_cw_or2_m;      dvs(cw_or2,:)];    
    
    dvs_cw_or3_r = [dvs_cw_or3_r; dvs_real(cw_or3,:)];    
    dvs_cw_or3_m = [dvs_cw_or3_m;      dvs(cw_or3,:)];    
    
    dvs_cw_or0_r = [dvs_cw_or0_r; dvs_real(cw_or0,:)];    
    dvs_cw_or0_m = [dvs_cw_or0_m;      dvs(cw_or0,:)];   
    
    % => ccw
    dvs_ccw_or1_r = [dvs_ccw_or1_r; dvs_real(ccw_or1,:)];    
    dvs_ccw_or1_m = [dvs_ccw_or1_m;      dvs(ccw_or1,:)];   
    
    dvs_ccw_or2_r = [dvs_ccw_or2_r; dvs_real(ccw_or2,:)];    
    dvs_ccw_or2_m = [dvs_ccw_or2_m;      dvs(ccw_or2,:)];   
    
    dvs_ccw_or3_r = [dvs_ccw_or3_r; dvs_real(ccw_or3,:)];    
    dvs_ccw_or3_m = [dvs_ccw_or3_m;      dvs(ccw_or3,:)];        
    
    % ==> index for computing correlations between higher res model fit and
    % raw DVs. This indexes the values in the model fit at the time points 
    % that correspond to the raw DV measures
    idx   = round(linspace(1,size(dvs_cw_or1_m,2),length(sacBegin:sacEnd))); %1:(size(dvs_cw_or1_m,2)/length(timeSac)):size(dvs_cw_or1_m,2);
    idx_r = sacBegin:sacEnd;
%     [val,pos]=intersect(timeSac,t);
    
    % ==> average over all individual DVs for each orientation (and for each session iSfln)
    dvs_cw_or1_m_mu = mean(dvs(cw_or1,:),1);
    dvs_cw_or1_r_mu = mean(dvs_real(cw_or1,:),1);
    orcw1r = corr(dvs_cw_or1_m_mu(idx)', dvs_cw_or1_r_mu(idx_r)','type','spearman');
    
    dvs_cw_or2_m_mu = mean(dvs(cw_or2,:),1);
    dvs_cw_or2_r_mu = mean(dvs_real(cw_or2,:),1);
    orcw2r = corr(dvs_cw_or2_m_mu(idx)', dvs_cw_or2_r_mu(idx_r)','type','spearman');   
    
    dvs_cw_or3_m_mu = mean(dvs(cw_or3,:),1);
    dvs_cw_or3_r_mu = mean(dvs_real(cw_or3,:),1);
    orcw3r = corr(dvs_cw_or3_m_mu(idx)', dvs_cw_or3_r_mu(idx_r)','type','spearman');     
    
    % ==> neutral vertical
    dvs_cw_or0_m_mu = mean(dvs(cw_or0,:),1);
    dvs_cw_or0_r_mu = mean(dvs_real(cw_or0,:),1);
    or0r = corr(dvs_cw_or0_m_mu(idx)', dvs_cw_or0_r_mu(idx_r)','type','spearman');  
    
    % ==> ccw 
    dvs_ccw_or1_m_mu = mean(dvs(ccw_or1,:),1);
    dvs_ccw_or1_r_mu = mean(dvs_real(ccw_or1,:),1);
    orccw1r = corr(dvs_ccw_or1_m_mu(idx)', dvs_ccw_or1_r_mu(idx_r)','type','spearman');  
    
    dvs_ccw_or2_m_mu = mean(dvs(ccw_or2,:),1);
    dvs_ccw_or2_r_mu = mean(dvs_real(ccw_or2,:),1);
    orccw2r = corr(dvs_ccw_or2_m_mu(idx)', dvs_ccw_or2_r_mu(idx_r)','type','spearman');     
    
    dvs_ccw_or3_m_mu = mean(dvs(ccw_or3,:),1);
    dvs_ccw_or3_r_mu = mean(dvs_real(ccw_or3,:),1);
    orccw3r = corr(dvs_ccw_or3_m_mu(idx)', dvs_ccw_or3_r_mu(idx_r)','type','spearman');    
    
    rs(iSfln,:) = [orcw1r, orcw2r, orcw3r, or0r, orccw1r, orccw2r, orccw3r]; 
 
end

%% =======

timeSac = S.dec.timeSac;
% Set time boundaries
fitTimeSacBegin = -800;  % in ms, relative to saccade onset
fitTimeSacEnd   = -50;
[~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
[~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));

% => figure 2 panel B
figure(); set(gcf,'color','white'); 
subplot(1,2,2); hold on; hold all;
plot(timeSac, mean(dvs_cw_or0_r,1),'.-', 'linewidth', 1,'markersize', 6, 'color', [1,1,1]*0.8); 
plot(timeSac, mean(dvs_cw_or1_r,1),'.-', 'linewidth', 2, 'markersize', 8, 'color', [1,0,0]*0.6); 
plot(timeSac, mean(dvs_cw_or2_r,1),'.-', 'linewidth', 3, 'markersize', 10, 'color', [1,0,0]*0.8); 
plot(timeSac, mean(dvs_cw_or3_r,1),'.-', 'linewidth', 4, 'markersize', 12, 'color', [1,0.2,0.2]); 

plot(timeSac, mean(dvs_ccw_or1_r,1),'b.-', 'linewidth', 2, 'markersize', 8, 'color', [0,0,1]*0.6); 
plot(timeSac, mean(dvs_ccw_or2_r,1),'b.-', 'linewidth', 3, 'markersize', 10, 'color', [0,0,1]*0.8); 
plot(timeSac, mean(dvs_ccw_or3_r,1),'b.-', 'linewidth', 4, 'markersize', 12, 'color', [0.2,0.2,1]); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
% xlim([timeSac(sacBegin),timeSac(sacEnd)]);
% => horizontal line at zero
plot(timeSac,zeros(1,size(timeSac,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[-0.6,0.6],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[-0.6,0.6],'k:')
ylim([-0.6,0.6]);
xlabel('Time');
ylabel('Signed Categorical DV');
title('Categorical Decision Variables (DVs)');
% legend('context 1','context 2')

% ==> model
% => figure
subplot(1,2,1); hold on; hold all; 
plot(t, mean(dvs_cw_or0_m,1),'r-');   
% => dv fit values at actual dv time intervals

plot(t, mean(dvs_cw_or0_m,1),'.-', 'linewidth', 1,'markersize', 4, 'color', [1,1,1]*0.8); 
plot(t, mean(dvs_cw_or1_m,1),'.-', 'linewidth', 2, 'markersize', 4, 'color', [1,0,0]*0.6); 
plot(t, mean(dvs_cw_or2_m,1),'.-', 'linewidth', 3, 'markersize', 4, 'color', [1,0,0]*0.8); 
plot(t, mean(dvs_cw_or3_m,1),'.-', 'linewidth', 4, 'markersize', 4, 'color', [1,0.2,0.2]); 

plot(t, mean(dvs_ccw_or1_m,1),'b.-', 'linewidth', 2, 'markersize', 4, 'color', [0,0,1]*0.6); 
plot(t, mean(dvs_ccw_or2_m,1),'b.-', 'linewidth', 3, 'markersize', 4, 'color', [0,0,1]*0.8); 
plot(t, mean(dvs_ccw_or3_m,1),'b.-', 'linewidth', 4, 'markersize', 4, 'color', [0.2,0.2,1]); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
% => horizontal line at zero
plot(t, zeros(1,size(dvs,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[-0.6,0.6],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[-0.6,0.6],'k:')
% xlim([timeSac(sacBegin),timeSac(sacEnd)]);
ylim([-0.6,0.6]);
xlabel('Time');
ylabel('Signed Categorical DV');
title('Model Fits');
% legend('context 1','context 2')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% => figure 2 panel B
figure(); set(gcf,'color','white'); 
subplot(1,2,2); hold on; hold all;
plot(timeSac, std(dvs_cw_or0_r,1),'.-', 'linewidth', 1,'markersize', 6, 'color', [1,1,1]*0.8); 
plot(timeSac, std(dvs_cw_or1_r,1),'.-', 'linewidth', 2, 'markersize', 8, 'color', [1,0,0]*0.6); 
plot(timeSac, std(dvs_cw_or2_r,1),'.-', 'linewidth', 3, 'markersize', 10, 'color', [1,0,0]*0.8); 
plot(timeSac, std(dvs_cw_or3_r,1),'.-', 'linewidth', 4, 'markersize', 12, 'color', [1,0.2,0.2]); 

plot(timeSac, std(dvs_ccw_or1_r,1),'b.-', 'linewidth', 2, 'markersize', 8, 'color', [0,0,1]*0.6); 
plot(timeSac, std(dvs_ccw_or2_r,1),'b.-', 'linewidth', 3, 'markersize', 10, 'color', [0,0,1]*0.8); 
plot(timeSac, std(dvs_ccw_or3_r,1),'b.-', 'linewidth', 4, 'markersize', 12, 'color', [0.2,0.2,1]); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
% xlim([timeSac(sacBegin),timeSac(sacEnd)]);
% => horizontal line at zero
plot(timeSac,zeros(1,size(timeSac,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[0,0.8],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[0,0.8],'k:')
ylim([0,0.8]);
xlabel('Time');
ylabel('Signed Categorical DV');
title('Categorical Decision Variables (DVs)');
% legend('context 1','context 2')

% time intervals
t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));

% ==> model
% => figure
subplot(1,2,1); hold on; hold all; 
plot(t, std(dvs_cw_or0_m,1),'r-');   
% => dv fit values at actual dv time intervals

plot(t, std(dvs_cw_or0_m,1),'.-', 'linewidth', 1,'markersize', 4, 'color', [1,1,1]*0.8); 
plot(t, std(dvs_cw_or1_m,1),'.-', 'linewidth', 2, 'markersize', 4, 'color', [1,0,0]*0.6); 
plot(t, std(dvs_cw_or2_m,1),'.-', 'linewidth', 3, 'markersize', 4, 'color', [1,0,0]*0.8); 
plot(t, std(dvs_cw_or3_m,1),'.-', 'linewidth', 4, 'markersize', 4, 'color', [1,0.2,0.2]); 

plot(t, std(dvs_ccw_or1_m,1),'b.-', 'linewidth', 2, 'markersize', 4, 'color', [0,0,1]*0.6); 
plot(t, std(dvs_ccw_or2_m,1),'b.-', 'linewidth', 3, 'markersize', 4, 'color', [0,0,1]*0.8); 
plot(t, std(dvs_ccw_or3_m,1),'b.-', 'linewidth', 4, 'markersize', 4, 'color', [0.2,0.2,1]); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
% => horizontal line at zero
plot(t, zeros(1,size(dvs,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[0,0.8],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[0,0.8],'k:')
% xlim([timeSac(sacBegin),timeSac(sacEnd)]);
ylim([0,0.8]);
xlabel('Time');
ylabel('Signed Categorical DV');
title('Model Fits');
% legend('context 1','context 2')

%% histogram of correlations (trial averaged (for each session), for each orientation (column dim) )

close all;

xlabs = [-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1];

% divided by monkey

rsF = rs(1:13,:);
rsJ = rs(14:end,:);

% edgs = [-1,-0.8,-0.6,-0.2,0,0.2,0.4,0.6,0.8,1];
edgs = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];


% ==> histogram of correlation coefficients
figure(); set(gcf,'color','white');
xticks(xlabs); xticklabels(xlabs); xlim([-1,1]); hold on; hold all;
histogram(rsF(:),edgs,'facecolor',[0.9,0.7,0]); hold on; hold all;
% plot([0,0],[0,1000],'r:','linewidth',3)
title('Correlation coefficients - Cat DV fit and real F');
ylabel('Frequency'); 
xlabel('Correlation')


% ==> histogram of correlation coefficients
figure(); set(gcf,'color','white');
xticks(xlabs); xticklabels(xlabs); xlim([-1,1]); hold on; hold all;
histogram(rsJ(:),edgs,'facecolor',[0,.7,.7]); hold on; hold all;
% plot([0,0],[0,1000],'r:','linewidth',3)
title('Correlation coefficients - Cat DV fit and real JP');
ylabel('Frequency');
xlabel('Correlation')

fprintf(['median animal F =', num2str(median(rsF(:))),'...\n'])
fprintf(['median animal JP =', num2str(median(rsJ(:))),'...\n'])


%% Updated F2D-E

bnd = 1;

timeSac = S.dec.timeSac;
% Set time boundaries
fitTimeSacBegin = -800;  % in ms, relative to saccade onset
fitTimeSacEnd   = -50;
[~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
[~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));

% => figure 2 panel B
figure(); set(gcf,'color','white'); 
subplot(1,2,1); hold on; hold all;
plot(timeSac, mean(dvs_cw_or0_r,1),'o:', 'linewidth', 1,'markersize', 8, 'color', [1,1,1]*0.8); 
plot(timeSac, mean(dvs_cw_or1_r,1),'o:', 'linewidth', 1, 'markersize', 8, 'color', [1,0,0]*0.6); 
plot(timeSac, mean(dvs_cw_or2_r,1),'o:', 'linewidth', 1, 'markersize', 8, 'color', [1,0,0]*0.8); 
plot(timeSac, mean(dvs_cw_or3_r,1),'o:', 'linewidth', 1, 'markersize', 8, 'color', [1,0.2,0.2]); 

plot(timeSac, mean(dvs_ccw_or1_r,1),'o:', 'linewidth', 1, 'markersize', 8, 'color', [0,0,1]*0.6); 
plot(timeSac, mean(dvs_ccw_or2_r,1),'o:', 'linewidth', 1, 'markersize', 8, 'color', [0,0,1]*0.8); 
plot(timeSac, mean(dvs_ccw_or3_r,1),'o:', 'linewidth', 1, 'markersize', 8, 'color', [0.2,0.2,1]); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
xlim([timeSac(sacBegin),timeSac(sacEnd)]);
% => horizontal line at zero
plot(timeSac,zeros(1,size(timeSac,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[-bnd, bnd],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[-bnd, bnd],'k:')
ylim([-bnd, bnd]);
xlabel('Time');
ylabel('Signed Categorical DV');
title('Categorical Decision Variables (DVs)');
% legend('context 1','context 2')

% time intervals
t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));

% ==> model
% => figure
subplot(1,2,1); hold on; hold all; 
plot(t, mean(dvs_cw_or0_m,1),'r-');   
% => dv fit values at actual dv time intervals

plot(t, mean(dvs_cw_or0_m,1),'.-', 'linewidth', 1,'markersize', 4, 'color', [1,1,1]*0.8); 
plot(t, mean(dvs_cw_or1_m,1),'.-', 'linewidth', 2, 'markersize', 4, 'color', [1,0,0]*0.6); 
plot(t, mean(dvs_cw_or2_m,1),'.-', 'linewidth', 3, 'markersize', 4, 'color', [1,0,0]*0.8); 
plot(t, mean(dvs_cw_or3_m,1),'.-', 'linewidth', 4, 'markersize', 4, 'color', [1,0.2,0.2]); 

plot(t, mean(dvs_ccw_or1_m,1),'b.-', 'linewidth', 2, 'markersize', 4, 'color', [0,0,1]*0.6); 
plot(t, mean(dvs_ccw_or2_m,1),'b.-', 'linewidth', 3, 'markersize', 4, 'color', [0,0,1]*0.8); 
plot(t, mean(dvs_ccw_or3_m,1),'b.-', 'linewidth', 4, 'markersize', 4, 'color', [0.2,0.2,1]); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
% => horizontal line at zero
plot(t, zeros(1,size(dvs,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[-bnd, bnd],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[-bnd, bnd],'k:')
xlim([timeSac(sacBegin),timeSac(sacEnd)]);
ylim([-bnd, bnd]);
xlabel('Time');
ylabel('Signed Categorical DV');
title('Model Fits');
% legend('context 1','context 2')




% ======> SI deviations plot

% ==> all single trial DV fits
ts_m = [dvs_cw_or0_m;  dvs_cw_or1_m;  dvs_cw_or2_m;  dvs_cw_or3_m; ...
        dvs_ccw_or1_m; dvs_ccw_or2_m; dvs_ccw_or3_m];
% ==> all single trial DVs
ts_r = [dvs_cw_or0_r;  dvs_cw_or1_r;  dvs_cw_or2_r;  dvs_cw_or3_r; ...
        dvs_ccw_or1_r; dvs_ccw_or2_r; dvs_ccw_or3_r]; 

% ==> index for model (model fit trajectories have 200 points for 40 DV points)
% idx = 1:(size(ts_m,2)/length(timeSac)):size(ts_m,2);

% % ==> square root of squared deviations between DV fits and real DVs
% % for each of the 40 timepoints
% err = sqrt((ts_m(:,idx) - ts_r).^2);

% ==> difference between DV fits and real DVs
% for each of the 40 timepoints
% err = ts_m(:,idx) - ts_r;

%% ==> F2A just real DV examples

close all; clc;

idxi = S.exp.taskContext == 1 ...
     & S.exp.stimOriDeg == 0 ...
     & S.exp.stimContrast == 1 ...
     & S.beh.choiceCat == 1 ... 
     & S.exp.respContext == 1;

ids = find(idxi == 1);

rd1 = 19; 
rd2 = 37;

% rd1 = randi(length(ids));
% rd2 = randi(length(ids));

% => sanity check -> visualize some examples
ex = ids(rd1);
corr(ts_m(ex,idx)',ts_r(ex,idx_r)','type','spearman')
figure(); set(gcf,'color','white');
plot(t,      ts_m(ex,:),'r-'); hold on; 
plot(timeSac,ts_r(ex,:),'ro:')
hold on; hold all;
ex = ids(rd2);
corr(ts_m(ex,idx)',ts_r(ex,idx_r)','type','spearman')
plot(t,      ts_m(ex,:),'b-'); hold on; 
plot(timeSac,ts_r(ex,:),'bo:');
xlim([timeSac(sacBegin),timeSac(sacEnd)]);
xticks(timeSac)
hold on; hold all;


% %% ==> F2A show some DV examples
% close all;
% 
% % => sanity check -> visualize some examples
% ex = 15279; %randi(size(ts_m,1),1);
% figure(); set(gcf,'color','white');
% plot(ts_m(ex,idx),'r.-'); hold on; plot(ts_r(ex,idx_r),'r+:','markerfacecolor','r')
% hold on; hold all;
% ex = 34170;
% plot(ts_m(ex,idx),'b.-'); hold on; plot(ts_r(ex,idx_r),'b+:','markerfacecolor','b')
% hold on; hold all;
% ex = 30100;
% plot(ts_m(ex,idx),'g.-'); hold on; plot(ts_r(ex,idx_r),'g+:','markerfacecolor','g')
% ex = 1590;
% plot(ts_m(ex,idx),'k.-'); hold on; plot(ts_r(ex,idx_r),'k+:','markerfacecolor','k')
% ylabel('Individual trial Categorical DVs')
% title('DV fits and real DVs')

% %% ==> plot deviations - a lot of points to draw
% 
% figure(); set(gcf,'color','white');
% plot(repmat(idx,[size(err,1),1]), err,'k.','markersize',0.2);
% hold on; hold all;
% plot(repmat(idx,[size(err,1),1]), mean(err,1),'k-')
% ylabel('distance between model fit point and actual DV point')
% set(gca,'xticklabels',[])

%% ==> rank correlations for each trial (DV fit & DV real)

edgs = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];

rks = nan(size(ts_m,1),1);

for i = 1:size(ts_m,1)
    % ==> store correlation coefficient
    % rks(i) = corr(ts_m(i,idx)', ts_r(i,:)');
    rks(i) = corr(ts_m(i,idx)', ts_r(i,idx_r)','type','spearman');
end

% ==> wilcoxon signed rank test
[P,H] = signrank(rks);

if iSfln == 13
    clrRGB = [0.9,0.7,0];
else
    clrRGB = [0,0.7,0.7];
end
% ==> histogram of correlation coefficients
fg = figure(); set(gcf,'color','white');  
% fg.Position = [647 549 598 330];
xticks(xlabs); xticklabels(xlabs); xlim([-1,1]); hold on; hold all;
histogram(rks,edgs,'facecolor',clrRGB); hold on; hold all;
% plot([0,0],[0,1000],'r:','linewidth',3)
title('Correlation coefficients - Cat DV fit and real');
ylabel('Frequency'); xlim([-1,1]);
xlabel(['Correlation, wilcoxon signed test p = ',num2str(P)])

fprintf(['median =', num2str(median(rks(:))),'...\n'])


%% == same plot just signed

% ts_m_cw  = [dvs_cw_or0_m;  dvs_cw_or1_m;  dvs_cw_or2_m;  dvs_cw_or3_m];
% ts_m_ccw = [dvs_ccw_or1_m; dvs_ccw_or2_m; dvs_ccw_or3_m];
% % ==> all single trial DVs
% ts_r_cw  = [dvs_cw_or0_r;  dvs_cw_or1_r;  dvs_cw_or2_r;  dvs_cw_or3_r];
% ts_r_ccw = [dvs_ccw_or1_r; dvs_ccw_or2_r; dvs_ccw_or3_r]; 
% 
% % == square root of squared deviations between DV fits and real DVs
% % for each of the 40 timepoints
% err_cw = sqrt((ts_m_cw(:,idx) - ts_r_cw).^2);
% 
% err_ccw = sqrt((ts_m_ccw(:,idx) - ts_r_ccw).^2);
% err_ccw = -err_ccw;
% 
% figure(); set(gcf,'color','white');
% plot(repmat(idx,[size(err_cw,1),1]),  err_cw, 'r.','markersize',0.2);
% hold on; hold all;
% plot(repmat(idx,[size(err_ccw,1),1]), err_ccw,'b.','markersize',0.2);
% hold on; hold all;
% plot(repmat(idx,[size(err_cw,1),1]),  mean(err_cw,1), 'r-')
% plot(repmat(idx,[size(err_ccw,1),1]), mean(err_ccw,1),'b-')
% ylabel('distance between model fit point and actual DV point')
% xlabel('time increment')
% set(gca,'xticklabels',[])


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Dir DV onset time and Cat DV peak time

clear all; close all; clc;
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/goris_data/data/';

ys  = nan(29,1);
SEs = nan(29,1);

for iSfln =  1:29 

% ===> Categorical DVs
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

% ===> Motor DVs
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% params: matrix with the model parameters [trial x parameter]
params_dir = load([drc,'trial_Dir_DV_params_iS_',num2str(iSfln),'.mat']);
params_dir = params_dir.ps_dir;
fprintf('Finished loading matrix with the model parameters for iS = %d... \n',iSfln)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% dvs are model predicted DV trajectories [trial x time] (Category DVs)
dvs_dir = load([drc,'trial_Dir_DV_traj_iS_',num2str(iSfln),'.mat']);
dvs_dir = dvs_dir.dv_dir;
fprintf('Finished loading model predicted DV trajectories (Directional DVs) for iS = %d... \n',iSfln)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> mini sanity check
% figure(); plot(mean(abs(dvs_dir),1),'b-'); hold on; hold all; plot(mean(abs(dvs),1),'r-')

% ==> load session data
if iSfln > 9
    load([dataPath,'/dataSet_',num2str(iSfln),'.mat']);
elseif iSfln <= 9
    load([dataPath,'/dataSet_0',num2str(iSfln),'.mat']);
end
% % ==> choice
% cho = S.beh.choiceCat;
% % ==> context, contrast, orientation indicator variables
% ctx = S.exp.taskContext; ctr = S.exp.stimContrast; ori = S.exp.stimOriDeg;

% ==> get times corresponding to peaks using indices in pks
% => from DV curve fit code...
timeSac = S.dec.timeSac;
% Set time boundaries
fitTimeSacBegin = -800;  % in ms, relative to saccade onset
fitTimeSacEnd   = -50;
[~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
[~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));
% ===> time intervals
t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~ Cat DV (signed) peak ~~~~~~~~~~~~~~~~~~~~~~~ %
% % ==> get max (maximum is either the positive peak if it's highest, or abs 
% % of negative peak if it's lowest)
% mx = [max(dvs, [], 2), abs(min(dvs, [], 2))];
% mxv = max(mx, [], 2);
% % => signed
% mxv(mx(:,1) <= mx(:,2)) = -mxv(mx(:,1) <= mx(:,2));
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~ Dir DV (signed) onset time ~~~~~~~~~~~~~~~~~~~~~ %
% ==> store index of first change in directional DV (using diff)
ontv = nan(size(dvs_dir,1),1);
for i = 1:size(dvs_dir,1)
    try
        %ontv(i) = params_dir(i,9) + params_dir(i,4) - 1.65*params_dir(i,3) + 1.65*params_dir(i,8);
        ontv(i) = params_dir(i,9) + params_dir(i,4) - 1.65*params_dir(i,3); % corrected ****
        
        % ==> index of first non-zero derivative (indicating a change) ****               
        %ontv(i) = find(diff(dvs_dir(i,:)) ~= 0, 1); 
        % => sanity check - visualize index/time of largest change
        %ontv(i) = find(diff(dvs_dir(i,:)) == max(diff(dvs_dir(i,:))), 1);
    catch
        fprintf('trial %d diff error\n...',i)
    end
end


% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% 
% % ==> get indicator matrix of Cat DV peaks
% pks = (dvs == mxv);
% 
% % keyboard
% 
% % ==> handle cases of multiple identical peaks (plateau)
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% % ==> indices of cases where there are multiple identical peaks
% pksi = find(sum(pks,2) > 1);
% % => cases where there are multiple identical peaks (plateau)
% m = pks(pksi,:);
% % => return indices for each row of first peak value 
% [~,i] = max(m,[],2);
% % => for rows with multiple peaks, just keep first peak index
% for j = 1:length(pksi)
%     idc = zeros(1,size(dvs,2)); 
%     idc(i(j)) = 1;
%     pks(pksi(j),:) = idc;
% end
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ==> Get Categorical DV peak times
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% => Cat DV (signed) peak times
CatDVts = nan(size(dvs,1),1);
for i = 1:size(dvs,1)
    % ==> get time from t intervals using index for trial i
    %CatDVts(i) = t(pks(i,:));
    CatDVts(i) = min(0, params(i, 4) + 1.65*params(i, 3)); %1.65
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% % ==> Motor Directional DV rise start times
% DirDVts = nan(size(dvs,1),1);
% for i = 1:size(dvs,1)
%     if ~isnan(ontv(i))
%         DirDVts(i) = t(ontv(i))';
%     end
% end
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% % the difference between the onset time of the MOTOR DV and 
% % peak time of the CATEGORICAL DV. A positive value indicates that the motor dv started 
% % to rise after the categorical dv peaked, and a negative value indicates that the motor 
% % dv started its excursion before the categorical dv peaked
% 
% % ==> differences
% dt = DirDVts(~isnan(ontv)) - CatDVts(~isnan(ontv));
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% ==> differences ****
DirDVts = ontv;
dt = DirDVts - CatDVts;

% ==> mean and SE
y = mean(dt);
% => standard error
SE = std(dt)/sqrt(length(dt));

% ==> store for all sessions - store dts in cell array but these should be
% a matrix...
ys(iSfln)  = y;
SEs(iSfln) = SE;

fprintf('completed differences for all trials for session %d of %d\n...',iSfln,29)

% % ==> sanity check -- view individual examples chosen at random
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% for i = 1:100
%     % ==> get indices of (signed) peaks (time indices)
%     tidx = randi(size(dvs,1));
%     % => peak time index
% %     pkt_idx = find(dvs(tidx,:) == mxv(tidx));
%     % ==> sometimes there's a plateau with multiple identical peaks
%     % => what then? -> choose first (earliest time) index
% %     if length(pkt_idx) > 1
% %         pkt_idx = pkt_idx(1);
% %         % ==> message
% %         fprintf('warning: multiple identical peaks...\n')
% %     end
%     % ==> mini sanity check
%     figure(1); set(gcf,'color','white');
%     plot(t,dvs(tidx,:),'r-'); hold all; %plot(t(pkt_idx), mxv(tidx),'ro'); hold on; hold all;
%     bndsc = [min(min(dvs(tidx,:),dvs_dir(tidx,:))),max(max(dvs(tidx,:),dvs_dir(tidx,:)))];
%     plot(repmat(CatDVts(tidx),1,2), bndsc,'r:'); hold on; hold all;
%     plot(t,dvs_dir(tidx,:),'b-'); hold all; 
%     % plot(t(ontv(tidx)), dvs_dir(tidx,ontv(tidx)),'bo'); hold on; hold all; ****
%     bndsd = [min(min(dvs(tidx,:),dvs_dir(tidx,:))),max(max(dvs(tidx,:),dvs_dir(tidx,:)))];
%     plot(repmat(DirDVts(tidx),1,2), bndsd,'b:'); hold on; hold all;
%     plot(t(1:end-1),diff(dvs_dir(tidx,:)),'k-')
%     %legend('Categorical DV','Cat DV peak','Dir DV','Dir DV rise start','Approx Derivative of Directional DV')
%     title(['time difference = ',num2str(dt(tidx))]);
%     drawnow; pause(10); 
%     clf;
% end
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

end

% ==> model parameters are:
% 1) - Bias in favor of CW or CCW choice
% 2) - Controls dynamic range of cat dv
% 3) - Controls speed rise of cat dv
% 4) - Controls half rise time of cat dv
% 5) - Controls decay after peak of cat dv
% 6) - Bias in favor of LW or RW choice
% 7) - Controls dynamic range of dir dv
% 8) - Controls speed rise of dir dv
% 9) - Controls half rise time of dir dv 



%% ==> plot
figure(2); grid on; set(gcf,'color','white');
for iSfln = 1:29
    title(['Difference (Dir offset and Cat peak) for all sessions']);
    axis square; hold on; hold all;
    ylabel('Session'); 
    yticks([0:30]); ylim([0,30]);
    errorbar(ys(iSfln),iSfln,SEs(iSfln),'r.-','horizontal')
    drawnow;
end

%% ==> dvCatPerf

dvCatPerf = load([drc,'dvCatPerf.mat'],'dvCatPerf');
dvCatPerf = dvCatPerf.dvCatPerf;

figure(); set(gcf,'color','white');
plot(dvCatPerf,ys,'o'); axis square;
xlabel('dvCatPerf');
ylabel('mean difference Dir DV onset time and Cat DV peak time')
title([ 'r = ', num2str(corr(dvCatPerf',ys,'type','spearman'))]);

text(dvCatPerf',ys,string(1:length(ys)))

[rs,ps] = corr(dvCatPerf',ys,'type','pearson')
[rm,pm] = corr(dvCatPerf',ys,'type','spearman')
 
%% F2 (old)

clear all; close all; clc;
dataPath = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/goris_data/data/';

% => one way is to store the indexed DVs across sessions
dvs_c1hi_r = [];
dvs_c1lo_r = [];
dvs_c2hi_r = [];
dvs_c2lo_r = [];

% => DV model fits
dvs_c1hi_m = [];
dvs_c1lo_m = [];
dvs_c2hi_m = [];
dvs_c2lo_m = [];

figure();
% ==> what is the session
for iSfln = 1:29 
    
    alpha_param = 'dynamic_range';

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

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % => single trial alpha values
    av = load([drc,'trial_DV_traj_alphas_iS_',num2str(iSfln),'_',alpha_param,'.mat'],  'av'); 
    av = av.av;
    % => indices (4 columns: 1=>idx_ccw_lo, 2=>idx_ccw_hi, 3=>idx_cw_lo, 4=> idx_cw_hi)
    idxm = load([drc,'trial_medSplit_idxVec_',num2str(iSfln),'_',alpha_param,'.mat'],  'idxm'); 
    idxm = idxm.idxm;
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
    ctx = S.exp.taskContext; 
    ctr = S.exp.stimContrast; 
    ori = S.exp.stimOriDeg;

    % => contrast values (changes depending on JP of FN!)
    ctrv = unique(ctr);
    % => unique orientations 
    oriv = unique(ori);

    % =-> actual Cat DVs
    dvs_real = S.dec.dvCatAll;

    % ==> subplots
    subplot(5,6,iSfln); hold on; hold all;
    
    % ==> low contrast, context 1 and 2 indices
    c1lo = (ctx==1)  & (ctr==ctrv(1));
    c2lo = (ctx==-1) & (ctr==ctrv(1));
    % => high contrast, context 1 and 2 indices
    c1hi = (ctx==1)  & (ctr==ctrv(2));
    c2hi = (ctx==-1) & (ctr==ctrv(2));
    %  ==> plot lo contrast cases
    plot(1:size(dvs_real,2),mean(dvs_real(c1lo,:),1),'r-'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(c2lo,:),1),'b-'); hold on; hold all;
    % ==> plot high contrast cases
    plot(1:size(dvs_real,2),mean(dvs_real(c1hi,:),1),'r:'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(c2hi,:),1),'b:'); hold on; hold all;
    ylim([-1,1]);
    drawnow;
    
    % ==> for F1 plot 
    dvs_c1lo_r = [dvs_c1lo_r; dvs_real(c1lo,:)];
    dvs_c1hi_r = [dvs_c1hi_r; dvs_real(c1hi,:)];
    
    dvs_c2lo_r = [dvs_c2lo_r; dvs_real(c2lo,:)];
    dvs_c2hi_r = [dvs_c2hi_r; dvs_real(c2hi,:)];  
    
    dvs_c1lo_m = [dvs_c1lo_m; dvs(c1lo,:)];
    dvs_c1hi_m = [dvs_c1hi_m; dvs(c1hi,:)];
    
    dvs_c2lo_m = [dvs_c2lo_m; dvs(c2lo,:)];
    dvs_c2hi_m = [dvs_c2hi_m; dvs(c2hi,:)];      
end

timeSac = S.dec.timeSac;
% Set time boundaries
fitTimeSacBegin = -800;  % in ms, relative to saccade onset
fitTimeSacEnd   = -50;

[~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
[~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));

% => figure 2 panel B
figure(); set(gcf,'color','white'); 
subplot(1,2,2); hold on; hold all;
plot(timeSac,mean(dvs_c1lo_r,1),'r.-','markersize',20); 
plot(timeSac,mean(dvs_c2lo_r,1),'b.-','markersize',20);
% ==> plot high contrast cases
plot(timeSac,mean(dvs_c1hi_r,1),'r.-'); 
plot(timeSac,mean(dvs_c2hi_r,1),'b.-'); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
xlim([min(timeSac),max(timeSac)]);
% => horizontal line at zero
plot(timeSac,zeros(1,size(timeSac,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[-0.3,0.3],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[-0.3,0.3],'k:')
ylim([-0.3,0.3]);
xlabel('Time');
ylabel('Signed Categorical DV');
title('Categorical Decision Variables (DVs)');
legend('context 1','context 2')

% time intervals
t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));

% ==> model
% => figure
subplot(1,2,1); hold on; hold all;
plot(t, mean(dvs_c1lo_m,1),'r-');  
plot(t, mean(dvs_c2lo_m,1),'b-');   
% => dv fit values at actual dv time intervals
scatter(t, mean(dvs_c1lo_m,1),20,'ro','filled');
scatter(t, mean(dvs_c2lo_m,1),20,'bo','filled');

% ==> plot high contrast cases
plot(t, mean(dvs_c1hi_m,1),'r-'); 
plot(t, mean(dvs_c2hi_m,1),'b-');
xticks(timeSac);
xticklabels(num2cell(timeSac))

% => dv fit values at actual dv time intervals
scatter(t, mean(dvs_c1hi_m,1),2,'ro','filled');
scatter(t, mean(dvs_c2hi_m,1),2,'bo','filled');
% => horizontal line at zero
plot(t, zeros(1,size(dvs,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[-0.3,0.3],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[-0.3,0.3],'k:')

xlim([min(timeSac),max(timeSac)]);
ylim([-0.3,0.3]);
xlabel('Time');
ylabel('Signed Categorical DV');
title('Model Fits');
legend('context 1','context 2')



%% F2A

iSfln = 17

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
ctx = S.exp.taskContext; 
ctr = S.exp.stimContrast; 
ori = S.exp.stimOriDeg;

% => contrast values (changes depending on JP of FN!)
ctrv = unique(ctr);
% => unique orientations 
oriv = unique(ori);

% ==> actual Cat DVs
dvs_real = S.dec.dvCatAll;


timeSac = S.dec.timeSac;
% Set time boundaries
fitTimeSacBegin = -800;  % in ms, relative to saccade onset
fitTimeSacEnd   = -50;

[~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
[~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));

% time intervals
t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));

ex = 15279; %randi(size(ts_m,1),1);
figure(); plot(ts_m(ex,idx),'ro-'); hold on; plot(ts_r(ex,:),'ko-')


%% ==> simulation F1

% => clear close, clear console
clear all; close all; clc;
% => relevant paths
dataPath     = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data';
drc = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/goris_data/data/';
pfc_functions_dr = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_analysis/simulation/pfc_functions';
addpath(pfc_functions_dr);

% verbosity for plots
plotFlag = 1;

% ==> markers for plots
mcs = {'k-','r-'};

% ==> use posterior or likelihood for plot?
pf = {'posterior','likelihood'};

for j = 1:2
    % => computing plots for likelihood or posterior
    f = pf{j}; 
    % => plot marker
    mc = mcs{j};
    % ==> run simulation and plotting
    IBOsimv2();

    % ==> bias and accuracy as a function of uncertainty
    fig4 = figure(6); fig4.set('Position', [971 469 888 349]);
    set(fig4,'color','white'); 
    subplot(1,2,1); hold on; hold all;
    plot(puv,bsv,mc,'linewidth',6); hold on; hold all;
    title('Bias and perceptual uncertainty');
    xlabel('perceptual uncertainty \sigma');
    ylabel('\Delta bias');
    legend('\Delta bias (posterior)','\Delta bias (likelihood)','Location','NorthWest')
    subplot(1,2,2); hold on; hold all;
    plot(puv,accs,mc,'linewidth',6); hold on; hold all;
    title('Choice accuracy and perceptual uncertainty')
    xlabel('perceptual uncertainty \sigma');
    ylabel('Overall accuracy'); 
    legend('Accuracy (posterior)','Accuracy (likelihood)');
end