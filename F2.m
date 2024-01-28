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



%% ==> histograms

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
for iSfln = 1:29 %1:13
    
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

%% 

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


% ==> all single trial DV fits
ts_m = [dvs_cw_or0_m;  dvs_cw_or1_m;  dvs_cw_or2_m;  dvs_cw_or3_m; ...
        dvs_ccw_or1_m; dvs_ccw_or2_m; dvs_ccw_or3_m];
% ==> all single trial DVs
ts_r = [dvs_cw_or0_r;  dvs_cw_or1_r;  dvs_cw_or2_r;  dvs_cw_or3_r; ...
        dvs_ccw_or1_r; dvs_ccw_or2_r; dvs_ccw_or3_r]; 


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


