%% ==> histograms

clear all; close all; clc;
dataPath = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '../data/';

% =========================
% ==> Monkey sessions
mIDs = 1:13;
% =========================

% => one way is to store the indexed DVs across sessions
dvs_cw_or1_r = [];
dvs_cw_or2_r = [];
dvs_cw_or3_r = [];
dvs_cw_or0_r = [];
dvs_ccw_or1_r = [];
dvs_ccw_or2_r = [];
dvs_ccw_or3_r = [];

% => DV model fits
dvs_cw_or1_m = [];
dvs_cw_or2_m = [];
dvs_cw_or3_m = [];
dvs_cw_or0_m = [];
dvs_ccw_or1_m = [];
dvs_ccw_or2_m = [];
dvs_ccw_or3_m = [];

rs = nan(29,7);

figure(); set(gcf,'Color','w'); 
% ==> what is the session
for iSfln = mIDs
    
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
    ctx = S.exp.taskContext; 
    ctr = S.exp.stimContrast; 
    ori = S.exp.stimOriDeg;

    % => contrast values (changes depending on JP of FN!)
    ctrv = unique(ctr);
    % => unique orientations 
    oriv = unique(ori); assert(length(oriv) == 7);

    % =-> actual Cat DVs
    dvs_real = S.dec.dvCatAll;

    % ==> subplots
    subplot(5,6,iSfln); hold on; hold all;

    % ==> high contrast, all orientations
    % => 1.1 degrees (JP)
    cw_or1 = (ori == oriv(end-2))           & (ctr==ctrv(2));
    % => 2.2 degrees (JP)
    cw_or2 = (ori == oriv(end-1))           & (ctr==ctrv(2));
    % => 3.3 degrees (JP)
    cw_or3 = (ori == oriv(end))             & (ctr==ctrv(2));
    % => 0 degrees
    cw_or0 = (ori == oriv(end-3))           & (ctr==ctrv(2));
    
    % => 11.1 degrees (JP)    
    ccw_or1 = (ori == oriv(3))              & (ctr==ctrv(2));
    % => -2.2 degrees (JP)
    ccw_or2 = (ori == oriv(2))              & (ctr==ctrv(2));
    % => -3.3 degrees (JP)
    ccw_or3 = (ori == oriv(1))              & (ctr==ctrv(2));
            
    % ==> plot 
    plot(1:size(dvs_real,2),mean(dvs_real(cw_or1,:),1),'r-'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(cw_or2,:),1),'r-'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(cw_or3,:),1),'r:'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(cw_or0,:),1),'r:'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(ccw_or1,:),1),'b-'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(ccw_or2,:),1),'b-'); hold on; hold all;
    plot(1:size(dvs_real,2),mean(dvs_real(ccw_or3,:),1),'b:'); hold on; hold all;
    
    ylim([-1,1]);
    drawnow;
    
    % ==> for F1 plot 
    % => cw
    dvs_cw_or1_r = [dvs_cw_or1_r; dvs_real(cw_or1,:)];    
    dvs_cw_or1_m = [dvs_cw_or1_m;      dvs(cw_or1,:)];
    
    dvs_cw_or2_r = [dvs_cw_or2_r; dvs_real(cw_or2,:)];    
    dvs_cw_or2_m = [dvs_cw_or2_m;      dvs(cw_or2,:)];    
    
    dvs_cw_or3_r = [dvs_cw_or3_r; dvs_real(cw_or3,:)];    
    dvs_cw_or3_m = [dvs_cw_or3_m;      dvs(cw_or3,:)];    
    
    % => vertical
    dvs_cw_or0_r = [dvs_cw_or0_r; dvs_real(cw_or0,:)];    
    dvs_cw_or0_m = [dvs_cw_or0_m;      dvs(cw_or0,:)];   
    
    % => ccw
    dvs_ccw_or1_r = [dvs_ccw_or1_r; dvs_real(ccw_or1,:)];    
    dvs_ccw_or1_m = [dvs_ccw_or1_m;      dvs(ccw_or1,:)];   
    
    dvs_ccw_or2_r = [dvs_ccw_or2_r; dvs_real(ccw_or2,:)];    
    dvs_ccw_or2_m = [dvs_ccw_or2_m;      dvs(ccw_or2,:)];   
    
    dvs_ccw_or3_r = [dvs_ccw_or3_r; dvs_real(ccw_or3,:)];    
    dvs_ccw_or3_m = [dvs_ccw_or3_m;      dvs(ccw_or3,:)];        
    
    % ==> index for computing correlations between higher res model fit and
    % raw DVs. This indexes the values in the model fit at the time points 
    % that correspond to the raw DV measures
    idx   = round(linspace(1,size(dvs_cw_or1_m,2),length(sacBegin:sacEnd))); %1:(size(dvs_cw_or1_m,2)/length(timeSac)):size(dvs_cw_or1_m,2);
    idx_r = sacBegin:sacEnd;
    
    % ==> average over all individual DVs for each orientation (and for each session iSfln)
    dvs_cw_or1_m_mu = mean(dvs(cw_or1,:),1);
    dvs_cw_or1_r_mu = mean(dvs_real(cw_or1,:),1);
    orcw1r = corr(dvs_cw_or1_m_mu(idx)', dvs_cw_or1_r_mu(idx_r)','type','spearman');
    
    dvs_cw_or2_m_mu = mean(dvs(cw_or2,:),1);
    dvs_cw_or2_r_mu = mean(dvs_real(cw_or2,:),1);
    orcw2r = corr(dvs_cw_or2_m_mu(idx)', dvs_cw_or2_r_mu(idx_r)','type','spearman');   
    
    dvs_cw_or3_m_mu = mean(dvs(cw_or3,:),1);
    dvs_cw_or3_r_mu = mean(dvs_real(cw_or3,:),1);
    orcw3r = corr(dvs_cw_or3_m_mu(idx)', dvs_cw_or3_r_mu(idx_r)','type','spearman');     
    
    % ==> neutral vertical
    dvs_cw_or0_m_mu = mean(dvs(cw_or0,:),1);
    dvs_cw_or0_r_mu = mean(dvs_real(cw_or0,:),1);
    or0r = corr(dvs_cw_or0_m_mu(idx)', dvs_cw_or0_r_mu(idx_r)','type','spearman');  
    
    % ==> ccw 
    dvs_ccw_or1_m_mu = mean(dvs(ccw_or1,:),1);
    dvs_ccw_or1_r_mu = mean(dvs_real(ccw_or1,:),1);
    orccw1r = corr(dvs_ccw_or1_m_mu(idx)', dvs_ccw_or1_r_mu(idx_r)','type','spearman');  
    
    dvs_ccw_or2_m_mu = mean(dvs(ccw_or2,:),1);
    dvs_ccw_or2_r_mu = mean(dvs_real(ccw_or2,:),1);
    orccw2r = corr(dvs_ccw_or2_m_mu(idx)', dvs_ccw_or2_r_mu(idx_r)','type','spearman');     
    
    dvs_ccw_or3_m_mu = mean(dvs(ccw_or3,:),1);
    dvs_ccw_or3_r_mu = mean(dvs_real(ccw_or3,:),1);
    orccw3r = corr(dvs_ccw_or3_m_mu(idx)', dvs_ccw_or3_r_mu(idx_r)','type','spearman');    
    
    rs(iSfln,:) = [orcw1r, orcw2r, orcw3r, or0r, orccw1r, orccw2r, orccw3r]; 
 
end

%% 
close all; clc;

timeSac = S.dec.timeSac;
% Set time boundaries
fitTimeSacBegin = -800;  % in ms, relative to saccade onset
fitTimeSacEnd   = -50;
[~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
[~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));

% => plot means
fg = figure(); set(fg,'color','white'); fg.Position = [133 572 1723 390];
subplot(1,4,2); hold on; hold all;
plot(timeSac, mean(dvs_cw_or0_r,1),'.-', 'linewidth', 1,'markersize', 6, 'color', [1,1,1]*0.8); 
plot(timeSac, mean(dvs_cw_or1_r,1),'.-', 'linewidth', 2, 'markersize', 8, 'color', [1,0,0]*0.6); 
plot(timeSac, mean(dvs_cw_or2_r,1),'.-', 'linewidth', 3, 'markersize', 10, 'color', [1,0,0]*0.8); 
plot(timeSac, mean(dvs_cw_or3_r,1),'.-', 'linewidth', 4, 'markersize', 12, 'color', [1,0.2,0.2]); 
% => ccw
plot(timeSac, mean(dvs_ccw_or1_r,1),'b.-', 'linewidth', 2, 'markersize', 8, 'color', [0,0,1]*0.6); 
plot(timeSac, mean(dvs_ccw_or2_r,1),'b.-', 'linewidth', 3, 'markersize', 10, 'color', [0,0,1]*0.8); 
plot(timeSac, mean(dvs_ccw_or3_r,1),'b.-', 'linewidth', 4, 'markersize', 12, 'color', [0.2,0.2,1]); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
xlim([timeSac(sacBegin),timeSac(sacEnd)]);
% => horizontal line at zero
plot(timeSac,zeros(1,size(timeSac,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[-0.6,0.6],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[-0.6,0.6],'k:')
ylim([-0.6,0.6]);
xlabel('Time');
ylabel('Signed Categorical DV (mean)');
title('Categorical DVs (mean)');

% ==> model
subplot(1,4,1); hold on; hold all; 
plot(t, mean(dvs_cw_or0_m,1),'r-');   
% => dv fit values at actual dv time intervals
plot(t, mean(dvs_cw_or0_m,1),'.-', 'linewidth', 1,'markersize', 4, 'color', [1,1,1]*0.8); 
plot(t, mean(dvs_cw_or1_m,1),'.-', 'linewidth', 2, 'markersize', 4, 'color', [1,0,0]*0.6); 
plot(t, mean(dvs_cw_or2_m,1),'.-', 'linewidth', 3, 'markersize', 4, 'color', [1,0,0]*0.8); 
plot(t, mean(dvs_cw_or3_m,1),'.-', 'linewidth', 4, 'markersize', 4, 'color', [1,0.2,0.2]); 
% => ccw
plot(t, mean(dvs_ccw_or1_m,1),'b.-', 'linewidth', 2, 'markersize', 4, 'color', [0,0,1]*0.6); 
plot(t, mean(dvs_ccw_or2_m,1),'b.-', 'linewidth', 3, 'markersize', 4, 'color', [0,0,1]*0.8); 
plot(t, mean(dvs_ccw_or3_m,1),'b.-', 'linewidth', 4, 'markersize', 4, 'color', [0.2,0.2,1]); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
% => horizontal line at zero
plot(t, zeros(1,size(dvs,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[-0.6,0.6],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[-0.6,0.6],'k:')
xlim([timeSac(sacBegin),timeSac(sacEnd)]);
ylim([-0.6,0.6]);
xlabel('Time');
ylabel('Model fit DV (mean)');
title('Model Fits (mean)');

% ==> plot standard deviations
subplot(1,4,4); hold on; hold all;
plot(timeSac, std(dvs_cw_or0_r,1),'.-', 'linewidth', 1,'markersize', 6, 'color', [1,1,1]*0.8); 
plot(timeSac, std(dvs_cw_or1_r,1),'.-', 'linewidth', 2, 'markersize', 8, 'color', [1,0,0]*0.6); 
plot(timeSac, std(dvs_cw_or2_r,1),'.-', 'linewidth', 3, 'markersize', 10, 'color', [1,0,0]*0.8); 
plot(timeSac, std(dvs_cw_or3_r,1),'.-', 'linewidth', 4, 'markersize', 12, 'color', [1,0.2,0.2]); 
% => ccw
plot(timeSac, std(dvs_ccw_or1_r,1),'b.-', 'linewidth', 2, 'markersize', 8, 'color', [0,0,1]*0.6); 
plot(timeSac, std(dvs_ccw_or2_r,1),'b.-', 'linewidth', 3, 'markersize', 10, 'color', [0,0,1]*0.8); 
plot(timeSac, std(dvs_ccw_or3_r,1),'b.-', 'linewidth', 4, 'markersize', 12, 'color', [0.2,0.2,1]); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
xlim([timeSac(sacBegin),timeSac(sacEnd)]);
% => horizontal line at zero
plot(timeSac,zeros(1,size(timeSac,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[0,0.8],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[0,0.8],'k:')
ylim([0,0.8]);
xlabel('Time');
ylabel('Signed Categorical DV (SD)');
title('Categorical DVs (SD)');

% time intervals
t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));

% ==> model
subplot(1,4,3); hold on; hold all; 
plot(t, std(dvs_cw_or0_m,1),'r-');   
% => dv fit values at actual dv time intervals
plot(t, std(dvs_cw_or0_m,1),'.-', 'linewidth', 1,'markersize', 4, 'color', [1,1,1]*0.8); 
plot(t, std(dvs_cw_or1_m,1),'.-', 'linewidth', 2, 'markersize', 4, 'color', [1,0,0]*0.6); 
plot(t, std(dvs_cw_or2_m,1),'.-', 'linewidth', 3, 'markersize', 4, 'color', [1,0,0]*0.8); 
plot(t, std(dvs_cw_or3_m,1),'.-', 'linewidth', 4, 'markersize', 4, 'color', [1,0.2,0.2]); 
% => ccw
plot(t, std(dvs_ccw_or1_m,1),'b.-', 'linewidth', 2, 'markersize', 4, 'color', [0,0,1]*0.6); 
plot(t, std(dvs_ccw_or2_m,1),'b.-', 'linewidth', 3, 'markersize', 4, 'color', [0,0,1]*0.8); 
plot(t, std(dvs_ccw_or3_m,1),'b.-', 'linewidth', 4, 'markersize', 4, 'color', [0.2,0.2,1]); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
% => horizontal line at zero
plot(t, zeros(1,size(dvs,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[0,0.8],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[0,0.8],'k:')
xlim([timeSac(sacBegin),timeSac(sacEnd)]);
ylim([0,0.8]);
xlabel('Time');
ylabel('Model fit DV (SD)');
title('Model Fits (SD)');


%% histogram of correlations (trial averaged (for each session)

close all; clc;

% ==> x axis labels
xlabs = -1:0.2:1;
% ==> edges
edgs = -1:0.1:1;

% ==> index monkey data
rsM = rs(mIDs,:);

% ==> histogram of correlation coefficients
fg = figure(); set(fg,'color','white'); 
xticks(xlabs); xticklabels(xlabs); xlim([-1,1]); hold on; hold all;
histogram(rsM(:),edgs,'facecolor',[0.9,0.7,0]); hold on; hold all;
title('Categorical DV fit and real DV Monkey F');
ylabel('Frequency'); 
xlabel('Correlation')

% ==> console update
fprintf(['median animal F =', num2str(median(rsM(:))),'...\n'])

%% Updated F2D-E

bnd = 1;

timeSac = S.dec.timeSac;
% Set time boundaries
fitTimeSacBegin = -800;  % in ms, relative to saccade onset
fitTimeSacEnd   = -50;
[~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
[~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));

% => figure 2 panel B
figure(); set(gcf,'color','white'); 
hold on; hold all;
plot(timeSac, mean(dvs_cw_or0_r,1),'o:', 'linewidth', 1,'markersize', 8, 'color', [1,1,1]*0.8); 
plot(timeSac, mean(dvs_cw_or1_r,1),'o:', 'linewidth', 1, 'markersize', 8, 'color', [1,0,0]*0.6); 
plot(timeSac, mean(dvs_cw_or2_r,1),'o:', 'linewidth', 1, 'markersize', 8, 'color', [1,0,0]*0.8); 
plot(timeSac, mean(dvs_cw_or3_r,1),'o:', 'linewidth', 1, 'markersize', 8, 'color', [1,0.2,0.2]); 
% => ccw
plot(timeSac, mean(dvs_ccw_or1_r,1),'o:', 'linewidth', 1, 'markersize', 8, 'color', [0,0,1]*0.6); 
plot(timeSac, mean(dvs_ccw_or2_r,1),'o:', 'linewidth', 1, 'markersize', 8, 'color', [0,0,1]*0.8); 
plot(timeSac, mean(dvs_ccw_or3_r,1),'o:', 'linewidth', 1, 'markersize', 8, 'color', [0.2,0.2,1]); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
xlim([timeSac(sacBegin),timeSac(sacEnd)]);
% => horizontal line at zero
plot(timeSac,zeros(1,size(timeSac,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[-bnd, bnd],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[-bnd, bnd],'k:')
ylim([-bnd, bnd]);
xlabel('Time');
ylabel('Signed Categorical DV');
title('Categorical Decision Variables (DVs)');


% time intervals
t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));

% ==> model
% => figure
plot(t, mean(dvs_cw_or0_m,1),'r-');   
% => dv fit values at actual dv time intervals
plot(t, mean(dvs_cw_or0_m,1),'.-', 'linewidth', 1,'markersize', 4, 'color', [1,1,1]*0.8); 
plot(t, mean(dvs_cw_or1_m,1),'.-', 'linewidth', 2, 'markersize', 4, 'color', [1,0,0]*0.6); 
plot(t, mean(dvs_cw_or2_m,1),'.-', 'linewidth', 3, 'markersize', 4, 'color', [1,0,0]*0.8); 
plot(t, mean(dvs_cw_or3_m,1),'.-', 'linewidth', 4, 'markersize', 4, 'color', [1,0.2,0.2]); 
% => ccw
plot(t, mean(dvs_ccw_or1_m,1),'b.-', 'linewidth', 2, 'markersize', 4, 'color', [0,0,1]*0.6); 
plot(t, mean(dvs_ccw_or2_m,1),'b.-', 'linewidth', 3, 'markersize', 4, 'color', [0,0,1]*0.8); 
plot(t, mean(dvs_ccw_or3_m,1),'b.-', 'linewidth', 4, 'markersize', 4, 'color', [0.2,0.2,1]); 
xticks(timeSac);
xticklabels(num2cell(timeSac))
% => horizontal line at zero
plot(t, zeros(1,size(dvs,2)),'k:')
plot([timeSac(sacBegin),timeSac(sacBegin)],[-bnd, bnd],'k:')
plot([timeSac(sacEnd),timeSac(sacEnd)],[-bnd, bnd],'k:')
xlim([timeSac(sacBegin),timeSac(sacEnd)]);
ylim([-bnd, bnd]);
xlabel('Time');
ylabel('Signed Categorical DV');
title('Model Fits');

%% ==> F2A just real DV examples

close all; clc;

% ==> all single trial DV fits
ts_m = [dvs_cw_or0_m;  dvs_cw_or1_m;  dvs_cw_or2_m;  dvs_cw_or3_m; ...
        dvs_ccw_or1_m; dvs_ccw_or2_m; dvs_ccw_or3_m];
% ==> all single trial DVs
ts_r = [dvs_cw_or0_r;  dvs_cw_or1_r;  dvs_cw_or2_r;  dvs_cw_or3_r; ...
        dvs_ccw_or1_r; dvs_ccw_or2_r; dvs_ccw_or3_r]; 


% ==> identify 'identical' trials.
idxi = S.exp.taskContext == 1 ...
     & S.exp.stimOriDeg == 0 ...
     & S.exp.stimContrast == 1 ...
     & S.beh.choiceCat == 1 ... 
     & S.exp.respContext == 1;

ids = find(idxi == 1);

% ==> nice curves for session 13 (monkey F)
rd1 = 19; 
rd2 = 37;

% ==> random choice
% rd1 = randi(length(ids));
% rd2 = randi(length(ids));

% => sanity check -> visualize some examples
ex = ids(rd1);
corr(ts_m(ex,idx)',ts_r(ex,idx_r)','type','spearman')
figure(); set(gcf,'color','white');
plot(t,      ts_m(ex,:),'r-'); hold on; 
plot(timeSac,ts_r(ex,:),'ro:')
hold on; hold all;
ex = ids(rd2);
corr(ts_m(ex,idx)',ts_r(ex,idx_r)','type','spearman')
plot(t,      ts_m(ex,:),'b-'); hold on; 
plot(timeSac,ts_r(ex,:),'bo:');
xlim([timeSac(sacBegin),timeSac(sacEnd)]);
xticks(timeSac)
hold on; hold all;


%% ==> rank correlations for each trial (DV fit & DV real)

% ==> array containing correlations
rks = nan(size(ts_m,1),1);

for i = 1:size(ts_m,1)
    % ==> store correlation coefficient
    rks(i) = corr(ts_m(i,idx)', ts_r(i,idx_r)','type','spearman');
end

% ==> wilcoxon signed rank test
[P,H] = signrank(rks);

% ==> plot colors
if iSfln == 13
    clrRGB = [0.9,0.7,0];
else
    clrRGB = [0,0.7,0.7];
end
% ==> histogram of correlation coefficients
fg = figure(); set(gcf,'color','white');  
% fg.Position = [647 549 598 330];
xticks(xlabs); xticklabels(xlabs); xlim([-1,1]); hold on; hold all;
histogram(rks,edgs,'facecolor',clrRGB); hold on; hold all;
title('Correlation coefficients - Cat DV fit and real');
ylabel('Frequency'); xlim([-1,1]);
xlabel(['Correlation, wilcoxon signed test p = ',num2str(P)])

% ==> print result
fprintf(['median = ', num2str(median(rks(:))),'...\n'])