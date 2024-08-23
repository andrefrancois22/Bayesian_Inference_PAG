clear all; close all; clc;
% rng(1,"twister");
% ==> Change to a local directory
FIGDR = 'Figures/';

% ==> parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> timepoints (like the one for Monkey F and session 1 in the physiology data)
t = linspace(-795,-45,200);
% get indices of timepoints closest to -500 and -300 ms
[~,il] = min(abs(t + 500)); %-500
[~,iu] = min(abs(t + 300)); %-300
% get indices of timepoints closest to -800 and -600 ms
[~,il2] = min(abs(t + 800)); %-800
[~,iu2] = min(abs(t + 600)); %-600

% => number of simulated trials
N = 3000;
tm = 200; % timepoints (same resolution as actual DV fits)  
% ==> bound
bnd = 8; %8

% ==> initial offset
ofs = 0.6; %0.6; 

% ==> accumulation bound or no bound?
M_FLAG = 'BOUND'; 
% M_FLAG = 'NO_BOUND';

% ==> add case 'IMPULSE_PRIOR'
% P_FLAG = 'IMPULSE_PRIOR';
P_FLAG = 'REG_DRIFT_PRIOR';

% ==> Trial-by-trial noise? ('cross-trial noise in the prior expectation')
% N_FLAG = 'NO_TRIAL_NOISE';
N_FLAG = 'TRIAL_NOISE';
% ==> define log normal noise
m = 1; % mean
v = 1; % variance

% ==> check

% ==> create a new SI that shows how the DDMs work:
% ==> show the prior step function and the sensory evidence 
% ==> also illustrate the cross-trial variability (the TRIAL_NOISE case)
% ==> where you show the multiplicative constant applied to the step
% function

% ==> to do
% have option to add noise to prior offset that varies across trials.

% ==> implementing contrast
% ==> try (1) change sd parameter 
% ==> try (2) change (reduce) intv

% ==> UPDATE METHODS - BASED ON IMPROVEMENTS IN ANALYSIS.
% ==> UPDATE FIGURE PANELS! ADD INTO THE OVERLEAF FOLDER! (F5 - x-axis values range should be a closer match)
% ==> F4A - add contrast manipulation PF curves.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> standard deviation for randn
sd = 1.25;
% ==> mean for randn
intv = 0.05;
mus = [-3*intv, -2*intv, -intv, 0, intv, 2*intv, 3*intv];
% mus = [-0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15];

% ==> prior offset factors (will influence drift rate due to integration)
fcs = 0.0:0.02:0.4; %0:0.01:0.2; 

% ==> run accumulation of evidence to bound drift diffusion (forward) model
% => compute simulated dynamic range split, and PFs
[predPF_ddm, db, dp, prop_cw, prop_ccw, dvs_c_cw, dvs_i_cw, dvs_c_ccw, dvs_i_ccw, csensb_cw_ds, csensb_ccw_ds] = accevbnd(N, tm, bnd, sd, mus, fcs, ofs, M_FLAG, P_FLAG, N_FLAG, m, v);

% ==> plot props or plot pfs?
plot_data = 'plot_props'; 
% plot_data = 'plot_pfs';

% ==> plot simulated choice proportions by dynamic range split
for fci = 3:length(fcs)
    fig3 = figure(3); set(gcf,'color','white'); set(gcf, 'Position',[273 218 1150 719]);
    subplot(4,7,fci)
    hold on; hold all; 
    plot(mus,squeeze(prop_cw(fci,3,:))', 'b-', 'linewidth', 1.5); 
    plot(mus,squeeze(prop_ccw(fci,3,:))','r-', 'linewidth', 1.5); 
    xlabel('Orientation');
    ylabel('p(cw)')    
    %xlim([-max(mus),max(mus)]); 
    ylim([0,1]);
    axis square;
    drawnow;      
    % ==> plot simulated choice proportions by dynamic range split
    fig4 = figure(4); set(gcf,'color','white'); set(gcf, 'Position',[273 218 1150 719]);
    subplot(4,7,fci)
    hold on; hold all; 
    if strcmp(plot_data,'plot_props')
        plot(mus,squeeze(prop_cw(fci,1,:))', 'color', [0.5,0.5,1], 'linewidth', 1.5); 
        plot(mus,squeeze(prop_ccw(fci,1,:))','color', [1,0.5,0.5],'linewidth', 1.5);  
        plot(mus,squeeze(prop_cw(fci,2,:))', 'b-', 'linewidth', 1.5); 
        plot(mus,squeeze(prop_ccw(fci,2,:))','r-', 'linewidth', 1.5); 
    end
    % ==> plot PFs over these
    if strcmp(plot_data,'plot_pfs')
        plot(mus,predPF_ddm{fci,1}(1,:),'b-','LineWidth',1.5)
        plot(mus,predPF_ddm{fci,1}(2,:),'r-','LineWidth',1.5)
        plot(mus,predPF_ddm{fci,2}(1,:),'color', [0.5,0.5,1],'LineWidth',1.5)
        plot(mus,predPF_ddm{fci,2}(2,:),'color', [1,0.5,0.5],'LineWidth',1.5)   
    end
    xlabel('Orientation');
    ylabel('p(cw)')
    %xlim([-max(mus),max(mus)]); 
    ylim([0,1]);
    axis square;
    drawnow;  
end
% ==> print figure to svg
print(fig4,[FIGDR,'SI-DYN_RANGE_PROPS', '-', P_FLAG, '-', M_FLAG, '-', N_FLAG, '-SD-', num2str(sd),'.svg'],'-dsvg')
% ==> print figure to svg
print(fig3,[FIGDR,'SI-SIMPLE_PROPS', '-', P_FLAG, '-', M_FLAG, '-', N_FLAG, '-SD-', num2str(sd),'.svg'],'-dsvg')

% ==> nice sanity check - plot all curves in same subpanels
% ==> plot PFs over these
for fci = 3:length(fcs)
    fig5 = figure(5); set(gcf,'color','white'); set(gcf, 'Position',[273 218 1150 719]);
    subplot(4,7,fci)
    hold on; hold all; 
    % plot(mus,squeeze(prop_cw(fci,3,:))', 'b-', 'linewidth', 1.5); 
    % plot(mus,squeeze(prop_ccw(fci,3,:))','r-', 'linewidth', 1.5); 
    % xlabel('Orientation');
    % ylabel('p(cw)')    
    %xlim([-max(mus),max(mus)]); 
    ylim([0,1]);
    axis square;
    plot(mus,squeeze(prop_cw(fci,1,:))', 'color', [0.5,0.5,1], 'linewidth', 1); 
    plot(mus,squeeze(prop_ccw(fci,1,:))','color', [1,0.5,0.5],'linewidth', 1);  
    plot(mus,squeeze(prop_cw(fci,2,:))', 'b-', 'linewidth', 1.5); 
    plot(mus,squeeze(prop_ccw(fci,2,:))','r-', 'linewidth', 1.5); 
    plot(mus,predPF_ddm{fci,1}(1,:),'bo-','LineWidth',1.5)
    plot(mus,predPF_ddm{fci,1}(2,:),'ro-','LineWidth',1.5)
    plot(mus,predPF_ddm{fci,2}(1,:),'color', [0.5,0.5,1],'LineWidth',1.5,'marker','o')
    plot(mus,predPF_ddm{fci,2}(2,:),'color', [1,0.5,0.5],'LineWidth',1.5,'marker','o')           
    xlabel('Orientation');
    ylabel('p(cw)')
    %xlim([-max(mus),max(mus)]); 
    ylim([0,1]);
    axis square;
    drawnow;  
end
% ==> print figure to svg
print(fig5,[FIGDR,'SI-PFS_AND_PROPS', '-', P_FLAG, '-', M_FLAG, '-', N_FLAG, '-SD-', num2str(sd),'.svg'],'-dsvg')


% ==> store averages of DVs at specific time windows
mu_early_csens_cw_ds = nan(length(fcs),1);
mu_late_csens_cw_ds  = nan(length(fcs),1);

mu_early_csens_ccw_ds = nan(length(fcs),1);
mu_late_csens_ccw_ds  = nan(length(fcs),1);

% ==> scatterplot for F3
for fci = 1:length(fcs)
    % ==> clockwise cases
    % ==> concatenate vectors across orientations
    csens_cw_dsc = vertcat(csensb_cw_ds{fci,:});
    mu_early_csens_cw_ds(fci)  = mean(mean(csens_cw_dsc(:,il2:iu2)));
    mu_late_csens_cw_ds(fci)   = mean(mean(csens_cw_dsc(:,il:iu)));
    % ==> counterclockwise cases
    % ==> concatenate vectors across orientations
    csens_ccw_dsc = vertcat(csensb_ccw_ds{fci,:});    
    mu_early_csens_ccw_ds(fci) = mean(mean(csens_ccw_dsc(:,il2:iu2)));
    mu_late_csens_ccw_ds(fci)  = mean(mean(csens_ccw_dsc(:,il:iu)));    
end

% ==> plot range
rg = floor(min([mu_early_csens_cw_ds; mu_early_csens_ccw_ds; mu_late_csens_cw_ds; mu_late_csens_ccw_ds])):ceil(max([mu_early_csens_cw_ds; mu_early_csens_ccw_ds; mu_late_csens_cw_ds; mu_late_csens_ccw_ds]));

figure(6); 
set(gcf,'color','white');
hold on; hold all;
scatter(mu_early_csens_cw_ds,mu_late_csens_cw_ds,   100, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'r');
scatter(mu_early_csens_ccw_ds,mu_late_csens_ccw_ds, 100, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'b');
plot(rg,rg,'k--');
plot(rg,zeros(length(rg)), 'k--');
plot(zeros(length(rg)),rg, 'k--');
xlabel('Early DV value (a.u)');
ylabel('Late DV value');
if strcmp(P_FLAG,'IMPULSE_PRIOR')
    title('Early bias + no drift');
elseif strcmp(P_FLAG, 'REG_DRIFT')
    title('Early bias + drift')
end
axis square;
drawnow;
% ==> print figure to svg
print(gcf,[FIGDR,'F3F-SCATTER', '-', P_FLAG, '-', M_FLAG, '-', N_FLAG, '-SD-', num2str(sd),'.svg'],'-dsvg')



% ==> correlation between simulated delta bias and simulated delta
% uncertainty
[r,p] = corr(db',dp');

% ==> plot db and dp
figure(7); set(gcf,'color','white');
scatter(db,dp,80,'ko','filled');
axis square;
xlabel('Simulated \Delta bias');
ylabel('Simulated \Delta perceptual uncertainty')
title(['r = ',num2str(r),', p = ',num2str(p)]);
% ==> print figure to svg
print(gcf,[FIGDR,'F5A-SCATTER', '-', P_FLAG, '-', M_FLAG, '-', N_FLAG, '-SD-', num2str(sd),'.svg'],'-dsvg')


% ==> Single congruent and incongruent DV
% ==> these may appear different from examples in the paper
ni = randi(size(dvs_i_ccw,1));
nc = randi(size(dvs_c_ccw,1)); 
% ==> draw figure
figure(8); set(gcf,'color','white');
hold on; hold all;
plot(t,dvs_i_ccw(ni,:)','linewidth',1.5, 'color', [0.6,0.6,0.6]);
plot(t,-bnd*ones(1,length(t)),'k--');
plot(t, bnd*ones(1,length(t)),'k--');
ylim([-(bnd+2), bnd+2]);
plot(t,dvs_c_ccw(nc,:)','linewidth',1.5,'color','k');
xlabel('Time');
ylabel('Accumulated evidence');
title('Example simulated trials');


%==> distribution of dynamic range trials
% ==> dynamic range function
dynf = @(x) max([max(x, [], 2) - x(:,1), -(min(x, [], 2) - x(:,1))], [], 2);

dynr_c_cw = dynf(dvs_c_cw);
dynr_i_cw = dynf(dvs_i_cw);

dynr_c_ccw = dynf(dvs_c_ccw); 
dynr_i_ccw = dynf(dvs_i_ccw); 

% ==> plot distribution of dynamic ranges
dynrs = [dynr_c_cw; dynr_i_cw; dynr_i_ccw; dynr_c_ccw];
figure(9); set(gcf,'color','white');
histogram(dynrs,15,'facecolor','w');
xlabel('Dynamic range (a.u)');
ylabel('Frequency');
title('Dynamic range distribution');
% ==> print figure to svg
print(gcf,[FIGDR,'F4C-DYN_RANGE_HIST', '-', P_FLAG, '-', M_FLAG, '-', N_FLAG, '-SD-', num2str(sd),'.svg'],'-dsvg')

