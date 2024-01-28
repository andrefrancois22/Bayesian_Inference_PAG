clear all; close all; clc;
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '../data/';

%% ==> compute DV peak averages by orientation

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

%% ==> Fit bias and slope to DV peaks

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

%% ==> Manuscript statistics for slope and bias differences

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

fg = figure(); set(fg,'color','white'); fg.Position = [139 604 1541 358];
subplot(1,4,1); 
scatter(slope(14:end,1),slope(14:end,2), 50,'o','filled', 'markeredgecolor', [1,0.5,0], 'markerfacecolor', [1,0.5,0]); axis square; ylim([0,0.6]); xlim([0,0.6]);
hold on; hold all;
plot(linspace(-0.1,0.65,10),linspace(-0.1,0.65,10),'k--')
xlabel('Slope low contrast')
ylabel('Slope high contrast')
title('JP')
subplot(1,4,2); 
scatter(slope(1:13,1),slope(1:13,2), 50,'o','filled', 'markeredgecolor', [0.15,0.75,0.5], 'markerfacecolor', [0.15,0.75,0.5]); axis square; ylim([0,0.4]); xlim([0,0.4]);
hold on; hold all;
plot(linspace(-0.1,0.65,10),linspace(-0.1,0.65,10),'k--')
xlabel('Slope low contrast')
ylabel('Slope high contrast')
title('F')

subplot(1,4,3); 
scatter(bv(14:end,1),bv(14:end,2), 50,'o','filled', 'markeredgecolor', [1,0.5,0], 'markerfacecolor', [1,0.5,0]); axis square; ylim([-0.4,1.1]); xlim([-0.4,1.1]);
hold on; hold all;
plot(linspace(-0.4,1,10),linspace(-0.4,1,10),'k--')
xlabel('Bias low contrast')
ylabel('Bias high contrast')
title('JP')
subplot(1,4,4); 
scatter(bv(1:13,1),bv(1:13,2), 50,'o','filled', 'markeredgecolor', [0.15,0.75,0.5], 'markerfacecolor', [0.15,0.75,0.5]); axis square; ylim([-0.1,0.65]); xlim([-0.1,0.65]);
hold on; hold all;
plot(linspace(-0.1,0.65,10),linspace(-0.1,0.65,10),'k--')
xlabel('Bias low contrast')
ylabel('Bias high contrast')
title('F')
