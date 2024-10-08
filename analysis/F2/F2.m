clear all; close all; clc;
% ==> directories
dataPath     = strcat('../../data/pfc_data/');
drc = '../../data/';

% ==> index vectors
idv  = cell(29,2,2,7);

% ==> raw DV values
dvsr = cell(29,2,2,7);
% ==> DV model fit values
dvsm = cell(29,2,2,7);

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
    timeSac = S.dec.timeSac;
    % Set time boundaries
    fitTimeSacBegin = -800;  % in ms, relative to saccade onset
    fitTimeSacEnd   = -50;
    [~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
    [~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));
    % time intervals
    t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));
    % ==> note: the t range is different for F and JP
    % => for F the range is   -795 -> -45
    % => for JP the range is  -790 -> -40
    
    % ==> context, contrast, orientation indicator variables
    ctx = S.exp.taskContext; ctr = S.exp.stimContrast; ori = S.exp.stimOriDeg;
    % x axis is orientation (the values differ for FN and JP!). Use unique values
    or = sort(unique(ori),'ascend')'; cx = sort(unique(ctx),'ascend')'; cr = sort(unique(ctr),'ascend')';            
    
    % ==> true dcCalAll (raw) values
    dvs_raw = S.dec.dvCatAll;
      
    for rc = 1:length(cr)
        for xc = 1:length(cx)
            for ro = 1:length(or)
                % ==> index vector
                idv{iS,xc,rc,ro} = (ori == or(ro)) & (ctx == cx(xc)) & (ctr == cr(rc));

                % ==> store raw DV values
                dvsr{iS,xc,rc,ro} = dvs_raw(idv{iS,xc,rc,ro},:);

                % ==> store model fit DV values
                dvsm{iS,xc,rc,ro} = dvs(idv{iS,xc,rc,ro},:);
            end
        end
    end
end

%%
clc;

% ==> orientation colormap
clrs = {[1.0,0.6,0.2],        ...
        [1.0,0.749,0.498],    ...
        [1.0,0.851,0.698],    ...
        [0.749,0.749,0.749],  ...
        [0.871,0.745,0.871],  ...
        [0.784,0.576,0.784], ...
        [0.655,0.325,0.655]};
% ==> draw figure
figure(); set(gcf,'color','white');

% ==> choose a session
for iS = 1:29

    % ==> load session data
    if iS > 9
        load([dataPath,'/dataSet_',num2str(iS),'.mat']);
    elseif iS <= 9
        load([dataPath,'/dataSet_0',num2str(iS),'.mat']);
    end
    timeSac = S.dec.timeSac;
    % Set time boundaries
    fitTimeSacBegin = -800;  % in ms, relative to saccade onset
    fitTimeSacEnd   = -50;
    [~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
    [~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));
    % time intervals
    t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));
    % ==> note: the t range is different for F and JP
    % => for F the range is   -795 -> -45
    % => for JP the range is  -790 -> -40
    
    % ==> orientation indicator variable
    ori = S.exp.stimOriDeg;
    % x axis is orientation (the values differ for FN and JP!). Use unique values
    or = sort(unique(ori),'ascend')';
    
    % ==> index for computing correlations between higher res model fit and
    % raw DVs. This indexes the values in the model fit at the time points 
    % that correspond to the raw DV measures
    idx   = round(linspace(1,size(t,2),length(sacBegin:sacEnd))); 
    idx_r = sacBegin:sacEnd;

    % ==> subplot for each session
    subplot(5,6,iS); 
    for ro = 1:length(or)
        % ==> model
        cmdv = vertcat(dvsm{iS,:,:,ro});
        % ==> raw
        crdv = vertcat(dvsr{iS,:,:,ro});
        % ==> plot
        hold on; hold all;
%         plot(t(idx),mean(cmdv(:,idx),1,'color', clrs{ro});
        plot(t,mean(cmdv,1),'color', clrs{ro});
        scatter(t(idx),mean(crdv(:,idx_r),1), 'markerfacecolor', clrs{ro},'markeredgecolor','w');
    end
    xlabel('Time from saccade (ms)');
    ylabel('Categorical DV');    
    % ==> title
    title([S.general.monkey,S.general.expDate])
    ylim([-1,1]);
    axis square;
    drawnow;
end

%%

% ==> plot averages for each monkey
M = {'F','J'};
rgs = {1:13,14:29};
% ==> draw figure (plotting average by monkey)
figure(); set(gcf,'color','white');

for m = 1:length(M)
    subplot(1,2,m)
    for ro = 1:length(or)
        % ==> model
        cmdv = vertcat(dvsm{rgs{m},:,:,ro});
        % ==> raw
        crdv = vertcat(dvsr{rgs{m},:,:,ro});
        % ==> plot
        hold on; hold all;
        
        fprintf('t must change here...\n')
        
        plot(t,mean(cmdv,1),'color', clrs{ro});
        scatter(t(idx),mean(crdv(:,idx_r),1), 'markerfacecolor', clrs{ro},'markeredgecolor','w');
    end
    xlabel('Time from saccade (ms)');
    ylabel('Categorical DV');
    % ==> title
    title(['Monkey ',M{m}])
    ylim([-1,1]);
    axis square;
    drawnow;    
end

%%

% ==> plot one example session
iS = 17; m = 2;

% ==> load session data
if iS > 9
    load([dataPath,'/dataSet_',num2str(iS),'.mat']);
elseif iS <= 9
    load([dataPath,'/dataSet_0',num2str(iS),'.mat']);
end
timeSac = S.dec.timeSac;
% Set time boundaries
fitTimeSacBegin = -800;  % in ms, relative to saccade onset
fitTimeSacEnd   = -50;
[~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
[~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));
% time intervals
t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));
% ==> note: the t range is different for F and JP
% => for F the range is   -795 -> -45
% => for JP the range is  -790 -> -40

% ==> orientation indicator variable
ori = S.exp.stimOriDeg;
% x axis is orientation (the values differ for FN and JP!). Use unique values
or = sort(unique(ori),'ascend')';

% ==> index for computing correlations between higher res model fit and
% raw DVs. This indexes the values in the model fit at the time points 
% that correspond to the raw DV measures
idx   = round(linspace(1,size(t,2),length(sacBegin:sacEnd))); 
idx_r = sacBegin:sacEnd;
    
% ==> draw figure (plotting average by monkey)
figure(); set(gcf,'color','white');

for ro = 1:length(or)
    % ==> model
    cmdv = vertcat(dvsm{iS,:,:,ro});
    % ==> raw
    crdv = vertcat(dvsr{iS,:,:,ro});
    % ==> plot
    hold on; hold all;
    plot(t,mean(cmdv,1),'color', clrs{ro});
    scatter(t(idx),mean(crdv(:,idx_r),1), 'markerfacecolor', clrs{ro},'markeredgecolor','w');
end
xlabel('Time from saccade (ms)');
ylabel('Categorical DV');
% ==> title
title([S.general.monkey,S.general.expDate])    
ylim([-1,1]);
axis square;
drawnow;  

%%
clc;

% ==> correlations
rs = cell(2,1);

% ==> correlation type
% rt = 'Spearman'; 
rt = 'Pearson';

% => Monkey colors
mclrs = {[1,0.75,0],[0.15,0.75,0.5]};

figure(); set(gcf,'color','white'); set(gcf,'Position',[738 380 1183 582]);
% ==> x axis labels
xlabs = -1:0.2:1;
% ==> edges
edgs = -1:0.1:1;

for m = 1:2
    for iS = rgs{m}

        % ==> load session data
        if iS > 9
            load([dataPath,'/dataSet_',num2str(iS),'.mat']);
        elseif iS <= 9
            load([dataPath,'/dataSet_0',num2str(iS),'.mat']);
        end
        timeSac = S.dec.timeSac;
        % Set time boundaries
        fitTimeSacBegin = -800;  % in ms, relative to saccade onset
        fitTimeSacEnd   = -50;
        [~, sacBegin] = min(abs(timeSac - fitTimeSacBegin));
        [~, sacEnd]   = min(abs(timeSac - fitTimeSacEnd));
        % time intervals
        t = linspace(timeSac(sacBegin),timeSac(sacEnd),size(dvs,2));
        % ==> note: the t range is different for F and JP
        % => for F the range is   -795 -> -45
        % => for JP the range is  -790 -> -40

        % ==> index for computing correlations between higher res model fit and
        % raw DVs. This indexes the values in the model fit at the time points 
        % that are closest to the raw DV measures (must round)
        idx   = round(linspace(1,size(t,2),length(sacBegin:sacEnd))); 
        idx_r = sacBegin:sacEnd;
        % ==> udpate
        fprintf('Finished loading session data for iS = %d... \n',iS)

        % ==> model
        cmdv = vertcat(dvsm{iS,:,:,:});
        % ==> raw
        crdv = vertcat(dvsr{iS,:,:,:});

        % ==> store correlations between raw DV values and corresponding DV
        % model values (matlab won't easily run corr on all trials for each animal, so append...)
        rs{m} = [rs{m}; diag(corr(cmdv(:,idx)',crdv(:,idx_r)','type',rt))];
    end
    subplot(1,2,m);
    xticks(xlabs); xticklabels(xlabs); xlim([-1,1]); hold on; hold all;
    histogram(rs{m},edgs,'facecolor',mclrs{m},'edgecolor','w'); hold on; hold all;    
    title(['Monkey ', M{m}, ' median ', rt,' \rho = ', num2str(median(rs{m}))]);
    xlabel([rt,' \rho']);
    ylabel('Frequency');
    axis square;
    drawnow;
end
