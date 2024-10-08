clear all; close all; clc;
% ==> directories
dataPath     = strcat('../../data/pfc_data/');
drc = '../../data/';

% ==> compute DV peak averages by orientation

% => initialize matrix containing values
mrx = nan(29,2,2,7);

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
    or = sort(unique(ori),'ascend')'; cx = sort(unique(ctx),'ascend')'; cr = sort(unique(ctr),'ascend')';
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % ~~~~~~~~~~~~~~~~~~~~~~~~ Cat DV (signed) peak ~~~~~~~~~~~~~~~~~~~~~~~ %
    % ==> get max (maximum is either the positive peak if it's highest, 
    % or abs of negative peak if it's lowest), but retain signed peaks.
    % **** NOTE: the initial offset value may be selected if its ****
    % **** absolute value is greater than any subsequent peaks   ****
    mx = [max(dvs, [], 2), abs(min(dvs, [], 2))];
    mxv = max(mx, [], 2);
    % => signed
    mxv(mx(:,1) <= mx(:,2)) = -mxv(mx(:,1) <= mx(:,2));
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    % => store averages for each stimulus condition
    % => loop through two context conditions
    for j = 1:length(cx)
        % => loop through two contrast conditions
        for k = 1:length(cr)
            % => loop through orientations
            for i = 1:length(or)
                mrx(iS,j,k,i) = mean( mxv( ori==or(i) & ctr==cr(k) & ctx==cx(j) ) );
            end                      
        end                    
    end      
    % => console updates...
    fprintf('Computed results for session %d of %d... \n',iS,29)
end

% ==> Fit bias and slope parameters to DV peaks
% => context colors
clrs = {[1,0,0],[0,0.75,1]};
    
% ==> bias and perceptual uncertainty vectors
bv  = nan(29,2);
puv = nan(29,2); 

% ==> for storing just slopes
slope = nan(29,2);
% => track loss for sanity check
lss = nan(29,2);

for iS =  1:29 
    
    % ==> draw figure
    figure(2); set(gcf,'color','white'); 
    clf; axis square;

    % => k is contrast
    for k = 1:2    

        % ==> fit linear model with identical slopes and different intercepts
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % ==> ys and x ir orientations
        y = squeeze(mrx(iS,:,k,:));
        x = or;
        % ==> fit linear models with identical slopes and different intercepts to
        % => initial parameter settings
        p0 = [0,0,0];
        % ==> two lines - 
        yh1 = @(p) (p(1)*x + p(2));
        yh2 = @(p) (p(1)*x + p(2) + p(3));
        % ==> MSE objective for both lines
        objective = @(p) sum( (yh1(p) - y(1,:)).^2 + (yh2(p) - y(2,:)).^2 );
        % ==> fmincon...
        p_opt = fmincon(objective,p0);
        % => inspect loss
        lss(iS,k) = objective(p_opt);        

        % ==> the bias is the vertical offset (optimal parameter 3 p(3))
        bv(iS,k)  = p_opt(3);
        % ==> uncertainty is slope
        puv(iS,k) = 1 / p_opt(1); %**** it's the inverse! ==> shallower slope is ~larger~ likelihood width      
        % ==> store slopes
        slope(iS,k) = p_opt(1);

        % ==> console update
        fprintf('completed linear model fit for context %d and contrast %k...\n',j,k)           
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        % ==> for visibility to inspect curve fits
        subplot(1,2,k)
        for j = 1:length(cx)  
            if j==1
                % ==> plot line fit          
                plot(or,yh1(p_opt),'color',clrs{j},'linewidth',3);
            elseif j==2
                % ==> plot line fit  
                plot(or,yh2(p_opt),'color',clrs{j},'linewidth',3);   
            end
            hold on;  hold all;
           plot(or,squeeze(mrx(iS,j,k,:))','o','markerfacecolor',[clrs{j}],'markeredgecolor','w','markersize',10); hold on; hold all;               
        end
        title('Peak Average (signed) by stimulus condition and orientation')    
        title(['Session: ',num2str(iS)]);
        ylabel('Peak Average (signed)');
        xlabel('Orientation');    
        ylim([-2,2]);    
    end
end

%%

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

% ==> Manuscript statistics for slope and bias differences

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


% => Monkey colors
mclrs = {[1,0.75,0],[0.15,0.75,0.5]};
% => monkey session ranges
rgs = {1:13,14:29};
% => monkeys
M = {'F','J'};

% ==> F2 - slope Hi v slope Lo, and bias Hi v bias Lo
fg = figure(); set(fg,'color','white'); fg.Position = [134 547 893 403];

for m = 1:length(M)
    subplot(1,2,1); 
    % ==> plot slope
    scatter(slope(rgs{m},1),slope(rgs{m},2), 150,'o','filled', 'markeredgecolor', 'w', 'markerfacecolor', mclrs{m}); axis square; %ylim([0,0.6]); xlim([0,0.6]);
    hold on; hold all;
    plot(linspace(-0.1,0.65,10),linspace(-0.1,0.65,10),'k--')
    plot(linspace(-0.1,0.65,10),zeros(1,length(linspace(-0.1,0.65,10))),'k--')
    plot(zeros(1,length(linspace(-0.1,0.65,10))),linspace(-0.1,0.65,10),'k--')      
    xlabel('Slope low contrast')
    ylabel('Slope high contrast')
    % ==> plot bias
    subplot(1,2,2); 
    scatter(bv(rgs{m},1),bv(rgs{m},2), 150,'o','filled', 'markeredgecolor', 'w', 'markerfacecolor', mclrs{m}); axis square; %ylim([-0.4,1.1]); xlim([-0.4,1.1]);
    hold on; hold all;
    plot(linspace(-0.4,1,10),linspace(-0.4,1,10),'k--')
    plot(linspace(-0.4,1,10),zeros(1,length(linspace(-0.4,1,10))),'k--')
    plot(zeros(1,length(linspace(-0.4,1,10))),linspace(-0.4,1,10),'k--')    
    xlabel('Bias low contrast')
    ylabel('Bias high contrast')
end

