clear all; close all; clc;
% ==> directories
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '../data/';

%% ==> accumulation bound panel

pr_cw  =  0.5;
pr_ccw = -0.5;

% => number of trials
N = 500;
% => number of timepoints
t = 200;

% ==>
sens  = randn(N,t);
% ==> cummulative
csens = cumsum(sens,2);

% ==> bound
bnd = 10;

% ==> csens with bound - two cases: with ccw prior or with cw prior
csensb_cw =  csens + pr_cw;
csensb_ccw = csens + pr_ccw;

% ==> set DV values for timepoints over trials that reach or exceed a bound
% to the bound value.
% ==> cw prior ctx
csensb_cw(csensb_cw >=  bnd) =  bnd;
csensb_cw(csensb_cw <= -bnd) = -bnd;
% ==> ccw prior ctx
csensb_ccw(csensb_ccw >=  bnd) =  bnd;
csensb_ccw(csensb_ccw <= -bnd) = -bnd;

% ==> once a bnd has been reached at time t for a given trial, set the remaining DV values
% starting from t+1 to the bound for that trial
% ==> cw cases
for i = 1:N
    bidx = min(find(csensb_cw(i,:) >=  bnd));
    % => set to bound
    if ~isempty(bidx)
        csensb_cw(i,bidx:end) = bnd;
    end
end
for i = 1:N
    bidx = min(find(csensb_cw(i,:) <=  -bnd));
    % => set to bound
    if ~isempty(bidx)
        csensb_cw(i,bidx:end) = -bnd;
    end
end
% ==> ccw context cases
for i = 1:N
    bidx = min(find(csensb_ccw(i,:) >=  bnd));
    % => set to bound
    if ~isempty(bidx)
        csensb_ccw(i,bidx:end) = bnd;
    end
end
for i = 1:N
    bidx = min(find(csensb_ccw(i,:) <=  -bnd));
    % => set to bound
    if ~isempty(bidx)
        csensb_ccw(i,bidx:end) = -bnd;
    end
end

figure(1); set(gcf,'Color','w'); set(gcf,'Position',[124 717 1682 245]);
subplot(1,5,1); 
hold on; hold all;
plot((csens + pr_cw)', 'color', [1,0.5,0.5,0.5]);     ylim([-bnd - 10, bnd + 10]);
plot(1:t,repmat( bnd,t), 'r--', 'linewidth', 4);
plot(1:t,repmat(-bnd,t), 'b--', 'linewidth', 4);
title('cw prior ctx')
subplot(1,5,2); 
hold on; hold all;
plot((csens + pr_ccw)', 'color', [0.5,0.5,1,0.5]);     ylim([-bnd - 10, bnd + 10]);
plot(1:t,repmat( bnd,t), 'r--', 'linewidth', 4);
plot(1:t,repmat(-bnd,t), 'b--', 'linewidth', 4);
title('ccw prior ctx')
 
subplot(1,5,3); 
hold on; hold all;
plot(1:t,repmat( bnd,t), 'r--', 'linewidth', 4);
plot(1:t,repmat(-bnd,t), 'b--', 'linewidth', 4);
plot(csensb_cw', 'color', [1,0.5,0.5,0.5]);     ylim([-bnd - 10, bnd + 10]);
title('cw prior ctx')
subplot(1,5,4); 
hold on; hold all;
plot(1:t,repmat( bnd,t), 'r--', 'linewidth', 4);
plot(1:t,repmat(-bnd,t), 'b--', 'linewidth', 4);
plot(csensb_ccw', 'color', [0.5,0.5,1,0.5]);     ylim([-bnd - 10, bnd + 10]);
title('ccw prior ctx')

subplot(1,5,5);
plot(mean(csensb_cw,1),'linewidth',1, 'color',  [1,0.5,0.5,0.5]);
hold on; hold all;
plot(mean(csensb_ccw,1),'linewidth',1, 'color', [0.5,0.5,1,0.5]);
ylabel('Average signed DV (acc. bound)')
xlabel('Time')
title('Average over trials')

%%

figure(); set(gcf,'Color','w');

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
    ctx = S.exp.taskContext; ctr = S.exp.stimContrast; ori = S.exp.stimOriDeg;
    
    % x axis is orientation (the values differ for FN and JP!). Use unique values
    or = unique(ori)'; cx = unique(ctx)'; cr = unique(ctr)';
    
    % indices
    idxcwlo =  (ori == 0) & (ctx == cx(2)) & (ctr == cr(1));
    idxccwlo = (ori == 0) & (ctx == cx(1)) & (ctr == cr(1));
    
    idxcwhi =  (ori == 0) & (ctx == cx(2)) & (ctr == cr(2));
    idxccwhi = (ori == 0) & (ctx == cx(1)) & (ctr == cr(2));    
    
    % ==> context cw (ctx == 1)
    vcwlo = dvs(idxcwlo,:);
    % ==> context cw (ctx == 1)
    vccwlo = dvs(idxccwlo,:);
    
    % ==> context cw (ctx == 1)
    vcwhi = dvs(idxcwhi,:);
    % ==> context cw (ctx == 1)
    vccwhi = dvs(idxccwhi,:);    
        
    
    subplot(5,6,iSfln);
    hold on; hold all;    
    plot(t,mean(vcwlo,1),  'r:');
    plot(t,mean(vccwlo,1), 'b:');
    plot(t,mean(vcwhi,1),  'r-');
    plot(t,mean(vccwhi,1), 'b-');    
    ylim([-1,1]);
    drawnow;
    
end

%% ==> just one example
close all; clc;

figure(); set(gcf,'Color','w');

iSfln = 17;

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
ctx = S.exp.taskContext; ctr = S.exp.stimContrast; ori = S.exp.stimOriDeg;

% x axis is orientation (the values differ for FN and JP!). Use unique values
or = unique(ori)'; cx = unique(ctx)'; cr = unique(ctr)';

% indices
idxcwlo =  (ori == 0) & (ctx == cx(2)) & (ctr == cr(1));
idxccwlo = (ori == 0) & (ctx == cx(1)) & (ctr == cr(1));

idxcwhi =  (ori == 0) & (ctx == cx(2)) & (ctr == cr(2));
idxccwhi = (ori == 0) & (ctx == cx(1)) & (ctr == cr(2));    

% ==> context cw (ctx == 1)
vcwlo = dvs(idxcwlo,:);
% ==> context cw (ctx == 1)
vccwlo = dvs(idxccwlo,:);

% ==> context cw (ctx == 1)
vcwhi = dvs(idxcwhi,:);
% ==> context cw (ctx == 1)
vccwhi = dvs(idxccwhi,:);    


subplot(1,2,1);
title(['J',S.general.expDate, ' low contrast']);
hold on; hold all;    
plot(t,mean(vcwlo,1),  'color', [1,0,0,0.5], 'linewidth',4);
plot(t,mean(vccwlo,1), 'color', [0,0,1,0.5], 'linewidth',4);
xlabel('Time'); ylabel('Average Signed DV (vertical stimulus)');
legend('CW prior context','CCW prior context');
ylim([-1,1]);
subplot(1,2,2);
title(['J',S.general.expDate, ' high contrast']);
hold on; hold all;   
plot(t,mean(vcwhi,1),  'r-','linewidth',4);
plot(t,mean(vccwhi,1), 'b-','linewidth',4);    
xlabel('Time'); ylabel('Average Signed DV (vertical stimulus)');
legend('CW prior context','CCW prior context');
ylim([-1,1]);
drawnow;
    
