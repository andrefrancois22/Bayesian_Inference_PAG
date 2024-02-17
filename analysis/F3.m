clear all; close all; clc;
% ==> directories
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '../data/';

pr_cw  =  0.5;
pr_ccw = -0.5;

% => number of trials
N = 1000;
% => number of timepoints
t = 100;

% ==>
sens  = randn(N,t);
% ==> cummulative
csens = cumsum(sens);
% ==> bound
bnd = 4;

% ==> csens with bound
csensb_cw = csens + pr_cw;

csensb_cw(csensb_cw >=  bnd) =  bnd;
csensb_cw(csensb_cw <= -bnd) = -bnd;

% figure(); 
% subplot(1,2,1); plot(csens,  'r');     ylim([-20,20]);
% subplot(1,2,2); plot(csensb_cw, 'r');  ylim([-20,20]);
% 
% figure();
% plot(mean(csensb_cw,1), 'r.-','linewidth',2);

%%

figure(); 

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
    plot(mean(vcwlo,1),  'r:');
    plot(mean(vccwlo,1), 'b:');
    plot(mean(vcwhi,1),  'r-');
    plot(mean(vccwhi,1), 'b-');    
    ylim([-1,1]);
    drawnow;
    
end