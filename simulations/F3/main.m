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
N = 10000;
tm = 200; % timepoints (same resolution as actual DV fits)  
% ==> bound
bnd = 12; 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> standard deviation for randn
sd = 1;
% ==> mean for randn
mus = [-0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15];% +0.08;
% ==> [mu_lb, mu_ub, intv, bs], (4 parameters) ==> linspace(mu_lb, mu_ub, intv) + bs

% ==> prior drift linear factors
fcs = 0:0.25:10; 
% (==> single parameter fc)

% ==> objective (NLL)
% ==> accevbnd outputs propCW for all 4 curves (from dynamic range split)
% iPF can range from 1-4, where propCW will contain all PFs (4 x 7 matrix of simulated proportions)
%  NLL(iPF) = -sum(log(max(1e-300, binopdf(nCW(iPF,:), nCCW(iPF,:) + nCW(iPF,:), propCW(iPF,:)))));

% ==> run accumulation of evidence to bound drift diffusion (forward) model
% => compute simulated dynamic range split, and PFs
[db,dp] = accevbnd(N, tm, bnd, sd, mus, fcs);


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

% ==> plot distribution of dynamic ranges
dynrs = [dynr_c_cw; dynr_i_cw; dynr_i_ccw; dynr_c_ccw];
figure(); set(gcf,'color','white');
histogram(dynrs,15,'facecolor','w');
xlabel('Dynamic range (a.u)');
ylabel('Frequency');
title('Dynamic range distribution');



%
% figure(); errorbar(1,mean(db),std(db),std(db))
% figure(); errorbar(1,mean(dp),std(dp),std(dp))