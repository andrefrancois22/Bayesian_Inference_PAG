

% ==> raw DV averages
cw_lo_raw_DVs  = [F_cw_lo_raw_DVs;  J_cw_lo_raw_DVs];
ccw_lo_raw_DVs = [F_ccw_lo_raw_DVs; J_ccw_lo_raw_DVs];

figure(); set(gcf,'Color','w'); set(gcf, 'Position',[675 363 602 599])

iSfln = 23;

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


title(['J',S.general.expDate]);
hold on; hold all;    
p1 = plot(t,mean(vcwlo,1),  'color', [1,0,0,0.5], 'linewidth',8);
p2 = plot(t,mean(vccwlo,1), 'color', [0,0,1,0.5], 'linewidth',8);
p4 = plot(timeSac(16:31), cw_lo_raw_DVs(iSfln,  16:31), 'o', 'color', [1,0,0,0.35]);
p5 = plot(timeSac(16:31), ccw_lo_raw_DVs(iSfln, 16:31), 'o', 'color', [0,0,1,0.35], 'linewidth',1);
% => time windows
xlabel('Time'); ylabel('Average Signed DV (vertical stimulus)');
ylim([-1,1]);
% => time windows
plot(repmat(-800,[20,1]),linspace(-0.8,0.8,20),'k--')
plot(repmat(-600,[20,1]),linspace(-0.8,0.8,20),'k--')
plot(repmat(-500,[20,1]),linspace(-0.8,0.8,20),'k:')
plot(repmat(-300,[20,1]),linspace(-0.8,0.8,20),'k:')
xlim([-840,-30]);

legend([p1,p2],'cw context, low contrast', ...
               'ccw context, low contrast');
drawnow;