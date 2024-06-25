clear all; close all; clc;
% ==> data directory
drc = '../../../data/';

% ==> dynamic range analysis types
dyn_types = {'popu','stim'};
% ==> load dynamic range results per stimulus orientation? or per population?
dyn_type = dyn_types{1};
% ==> load choice proportions and counts from dynamic range split on real data
% CPs (<Session ID>_<contrast>_<dynamic range split>)
load([drc,'CPs_dynr_',dyn_type,'.mat'],'CPs');
% CTs (<Session ID>_<contrast>_<dynamic range split>_<counts (nccw (1) and ncw (2))>)
load([drc,'CTs_dynr_',dyn_type,'.mat'],'CTs');
% ==> load PFs
% PFs (<Session ID>_<contrast>_<dynamic range split>)
load([drc,'PFs_dynr_',dyn_type,'.mat'],'PF');

% ==> parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% => number of simulated trials
N = 1000;
tm = 200; % timepoints (same resolution as actual DV fits)  
% ==> bound
bnd = 12; 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> standard deviation for randn
sd = 1;
% ==> mean for randn
mus = [-0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15] + 0.1;
% ==> [intv, bs], (2 parameters) ==> [-3*intv, -2*intv, -intv, 0, intv, 2*intv, 3*intv] + bs

% ==> prior drift linear factors
fc = 5; 
% (==> single parameter fc)
ofs = 15;

% ==> run accumulation of evidence to bound drift diffusion (forward) model
% => compute simulated dynamic range split, and PFs
[propcw, propccw] = accevbndf(N, tm, bnd, sd, mus, fc, ofs);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ==> PF curve estimation Set options
fitFlag         = 1;       % 0 or 1, where 1 means do the fit right now
nBootStraps     = 1000;     % determines number of non parametric bootstraps if fitFlag == 1
options         = optimoptions('fmincon','ConstraintTolerance', 1e-100, ...
                                         'OptimalityTolerance', 1e-100, ...
                                         'StepTolerance',       1e-100, ...
                                         'Algorithm', 'sqp');
options.Display = 'iter';
% Set bounds on model parameters
% Model fit to choice data split by task context (1 & 2: guess rate; 3: perceptual uncertainty; 4 & 5: decision criterion)
% ==> interval, overall bias, fc (drift rate), initial offset, bound
startVec = [0.038, 0, 8, 90];%, 12]; 
LB(1,1)  = 0.01;                      UB(1,1) = 0.06;    % orientation mean interval
LB(2,1)  = -0.2;                      UB(2,1) = 0.2;  % overall bias (not decision bias - this is left or right shift of all curves)
LB(3,1)  = 0;                         UB(3,1) = 20;   % drift rate (slope of the linear drift)
LB(4,1)  = 0;                         UB(4,1) = 200;   % initial offset
% LB(5,1)  = 6;                         UB(5,1) = 20;   % bound
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> for each session, and for each contrast level, fit acc-bound predicted dynr PFs curves to actual data
iS = 17;
cr = 2; % => lo contrast (fit models separately for different stimulus contrast images)
% ==> wrap in objective function
[nll] = calcnllf(iS,cr,CTs,propccw,propcw);
fprintf('NLL = %d...\n',nll)

mod_type = 'm2';

% ==> objective function for fmincon
obFun = @(paramVec) modfitf(iS, cr, CTs, N, tm, sd, paramVec, mod_type); 
% ==> run optimization
simfit = fmincon(obFun, startVec, [], [], [], [], LB, UB, [], options);
% ==> return best fit simulated dynamic range proportions
[nll, propcw, propccw] = modfitf(iS, cr, CTs, N, tm, sd, simfit, mod_type);

simfit
nll
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

% simfit = [0.038, 0, 8, 90]; mod_type = 'm2'; %(17 ctr 1)
simfit = [0.05, 0, 3, 40]; mod_type = 'm2';    %(17, ctr 2)
[nll, propcw, propccw] = modfitf(iS, cr, CTs, N, tm, sd, simfit, mod_type);
nll

close all; 
figure(1);
subplot(2,3,1);
hold on; hold all;
plot(1:7,propcw(1,:),'r-h')
plot(1:7,propccw(1,:),'b-h')
plot(1:7,propcw(2,:),'rh:')
plot(1:7,propccw(2,:),'bh:')
subplot(2,3,2);
hold on; hold all;
plot(CPs{iS,cr,2}(2,:),'b:')
plot(CPs{iS,cr,2}(1,:),'r:')
plot(CPs{iS,cr,1}(1,:),'r')
plot(CPs{iS,cr,1}(2,:),'b')
subplot(2,3,3);
hold on; hold all;
plot(1:7, CTs{iS,cr,1,2}(2,:) ./ (CTs{iS,cr,1,2}(2,:) + CTs{iS,cr,1,1}(2,:)),'b');
plot(1:7, CTs{iS,cr,1,2}(1,:) ./ (CTs{iS,cr,1,2}(1,:) + CTs{iS,cr,1,1}(1,:)),'r');
plot(1:7, CTs{iS,cr,2,2}(2,:) ./ (CTs{iS,cr,2,2}(2,:) + CTs{iS,cr,2,1}(2,:)),'b:');
plot(1:7, CTs{iS,cr,2,2}(1,:) ./ (CTs{iS,cr,2,2}(1,:) + CTs{iS,cr,2,1}(1,:)),'r:');
hold on; hold all;
plot(1:7,propcw(1,:), 'r-h')
plot(1:7,propccw(1,:),'b-h')
plot(1:7,propcw(2,:), 'rh:')
plot(1:7,propccw(2,:),'bh:')


subplot(2,3,5);
hold  on; hold all;
plot(1:7,propcw(2,:), 'rh:')
plot(1:7,propccw(2,:),'bh:')
plot(1:7,PF{iS,cr,2}(2,:),'b-')
plot(1:7,PF{iS,cr,2}(1,:),'r-')
subplot(2,3,6);
hold on; hold all;
plot(1:7,propcw(1,:), 'r-h')
plot(1:7,propccw(1,:),'b-h')
plot(1:7,PF{iS,cr,1}(2,:),'b-')
plot(1:7,PF{iS,cr,1}(1,:),'r-')


