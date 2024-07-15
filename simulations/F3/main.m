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
N = 1000;
tm = 200; % timepoints (same resolution as actual DV fits)  
% ==> bound
bnd = 12;
% ==> initial offset
ofs = 0.25; 

% M_FLAG = 'BOUND'; 
M_FLAG = 'NO_BOUND';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> standard deviation for randn
sd = 1;
% ==> mean for randn
mus = [-0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15];

% ==> prior offset factors (will influence drift rate due to integration)
fcs = 0:0.01:0.4; 

% ==> run accumulation of evidence to bound drift diffusion (forward) model
% => compute simulated dynamic range split, and PFs
[db,dp, prop_cw, prop_ccw, dvs_c_cw, dvs_i_cw, dvs_c_ccw, dvs_i_ccw] = accevbnd(N, tm, bnd, sd, mus, fcs, ofs, M_FLAG);

%%

% ==> plot simulated choice proportions by dynamic range split
for fci = 1:length(fcs)
    figure(3); set(gcf,'color','white'); set(gcf, 'Position',[273 218 1150 719]);
    subplot(6,7,fci)
    hold on; hold all; 
    plot(mus,squeeze(prop_cw(fci,3,:))','b.-', 'linewidth', 1.5); 
    plot(mus,squeeze(prop_ccw(fci,3,:))','r.-', 'linewidth', 1.5); 
    xlabel('Orientation');
    ylabel('p(cw)')    
    %xlim([-max(mus),max(mus)]); 
    ylim([0,1]);
    axis square;
    drawnow;      
    % ==> plot simulated choice proportions by dynamic range split
    figure(4); set(gcf,'color','white'); set(gcf, 'Position',[273 218 1150 719]);
    subplot(6,7,fci)
    hold on; hold all; 
    plot(mus,squeeze(prop_cw(fci,1,:))','bx:', 'linewidth', 1); 
    plot(mus,squeeze(prop_ccw(fci,1,:))','rx:','linewidth', 1);  
    plot(mus,squeeze(prop_cw(fci,2,:))','b.-', 'linewidth', 1.5); 
    plot(mus,squeeze(prop_ccw(fci,2,:))','r.-', 'linewidth', 1.5); 
    xlabel('Orientation');
    ylabel('p(cw)')
    %xlim([-max(mus),max(mus)]); 
    ylim([0,1]);
    axis square;
    drawnow;  
end

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
ni = 613; %randi(size(dvs_i_ccw,1));
nc = 91;  %randi(size(dvs_c_ccw,1)); 

% ==> draw figure
figure(); set(gcf,'color','white');
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
