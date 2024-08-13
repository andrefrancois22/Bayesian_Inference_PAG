clear all; close all; clc;

rng(1,"twister");

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
N = 2500;
tm = 200; % timepoints (same resolution as actual DV fits)  
% ==> bound
bnd = 8; %12; -6
% ==> initial offset
ofs = 0.6; %0.25; 

% ==> accumulation bound or no bound?
% M_FLAG = 'BOUND'; 
M_FLAG = 'NO_BOUND';

% ==> add case 'IMPULSE_PRIOR'
% P_FLAG = 'IMPULSE_PRIOR';
P_FLAG = 'REG_DRIFT_PRIOR';

% ==> Trial-by-trial noise? ('cross-trial noise in the prior expectation')
N_FLAG = 'NO_TRIAL_NOISE';
% N_FLAG = 'TRIAL_NOISE';

% ==> check

% ==> create a new SI that shows how the DDMs work:
% ==> show the prior step function and the sensory evidence 
% ==> also illustrate the cross-trial variability (the TRIAL_NOISE case)
% ==> where you show the multiplicative constant applied to the step
% function

% ==> to do
% have option to add noise to prior offset that varies across trials.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> standard deviation for randn
sd = 1;
% ==> mean for randn
intv = 0.05;
mus = [-3*intv, -2*intv, -intv, 0, intv, 2*intv, 3*intv];
% mus = [-0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15];

% ==> prior offset factors (will influence drift rate due to integration)
fcs = 0.0:0.02:0.4; %0:0.01:0.2; 

% ==> run accumulation of evidence to bound drift diffusion (forward) model
% => compute simulated dynamic range split, and PFs
[db,dp, prop_cw, prop_ccw, dvs_c_cw, dvs_i_cw, dvs_c_ccw, dvs_i_ccw, csensb_cw_ds, csensb_ccw_ds] = accevbnd(N, tm, bnd, sd, mus, fcs, ofs, M_FLAG, P_FLAG, N_FLAG);


% ==> plot simulated choice proportions by dynamic range split
for fci = 3:length(fcs)
    figure(3); set(gcf,'color','white'); set(gcf, 'Position',[273 218 1150 719]);
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
    figure(4); set(gcf,'color','white'); set(gcf, 'Position',[273 218 1150 719]);
    subplot(4,7,fci)
    hold on; hold all; 
    plot(mus,squeeze(prop_cw(fci,1,:))', 'color', [0.5,0.5,1], 'linewidth', 1.5); 
    plot(mus,squeeze(prop_ccw(fci,1,:))','color', [1,0.5,0.5],'linewidth', 1.5);  
    plot(mus,squeeze(prop_cw(fci,2,:))', 'b-', 'linewidth', 1.5); 
    plot(mus,squeeze(prop_ccw(fci,2,:))','r-', 'linewidth', 1.5); 
    xlabel('Orientation');
    ylabel('p(cw)')
    %xlim([-max(mus),max(mus)]); 
    ylim([0,1]);
    axis square;
    drawnow;  
end

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

figure(10); 
set(gcf,'color','white');
hold on; hold all;
scatter(mu_early_csens_cw_ds,mu_late_csens_cw_ds,   100, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'r');
scatter(mu_early_csens_ccw_ds,mu_late_csens_ccw_ds, 100, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'b');
plot(rg,rg,'k--');
plot(rg,zeros(length(rg)), 'k--');
plot(zeros(length(rg)),rg, 'k--');
xlabel('Early DV value (a.u)');
ylabel('Late DV value');
title('Early bias + drift');
axis square;
drawnow;

%% 
% ==> correlation between simulated delta bias and simulated delta
% uncertainty
[r,p] = corr(db',dp');

% ==> plot db and dp
figure(6); set(gcf,'color','white');
scatter(db,dp,80,'ko','filled');
axis square;
xlabel('Simulated \Delta bias');
ylabel('Simulated \Delta perceptual uncertainty')
title(['r = ',num2str(r),', p = ',num2str(p)]);

%% 
% ==> Single congruent and incongruent DV

% ==> these may appear different from examples in the paper
ni = randi(size(dvs_i_ccw,1));
nc = randi(size(dvs_c_ccw,1)); 

% ==> draw figure
figure(3); set(gcf,'color','white');
hold on; hold all;
plot(t,dvs_i_ccw(ni,:)','linewidth',1.5, 'color', [0.6,0.6,0.6]);
plot(t,-bnd*ones(1,length(t)),'k--');
plot(t, bnd*ones(1,length(t)),'k--');
ylim([-(bnd+2), bnd+2]);
plot(t,dvs_c_ccw(nc,:)','linewidth',1.5,'color','k');
xlabel('Time');
ylabel('Accumulated evidence');
title('Example simulated trials');

%% 
%==> distribution of dynamic range trials

% ==> dynamic range function
dynf = @(x) max([max(x, [], 2) - x(:,1), -(min(x, [], 2) - x(:,1))], [], 2);

dynr_c_cw = dynf(dvs_c_cw);
dynr_i_cw = dynf(dvs_i_cw);

dynr_c_ccw = dynf(dvs_c_ccw); 
dynr_i_ccw = dynf(dvs_i_ccw); 

% ==> plot distribution of dynamic ranges
dynrs = [dynr_c_cw; dynr_i_cw; dynr_i_ccw; dynr_c_ccw];
figure(); set(gcf,'color','white');
histogram(dynrs,15,'facecolor','w');
xlabel('Dynamic range (a.u)');
ylabel('Frequency');
title('Dynamic range distribution');
