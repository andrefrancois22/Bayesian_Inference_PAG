clear all; close all; clc;
% ==> directories
dataPath     = strcat('../../data/pfc_data/');
drc = '../../data/';

% ==> draw figure
figure(); set(gcf,'Color','w');

% ==> index vectors
idv  = cell(29,2,2);
% ==> dv values in early and late time windows
dvw  = cell(29,2,2);
% ==> raw DV values
dvsr = cell(29,2,2);
% ==> DV model fit values
dvsm = cell(29,2,2);

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
    
    % get indices of timepoints closest to -500 and -300 ms
    [~,il] = min(abs(t + 500)); %-500
    [~,iu] = min(abs(t + 300)); %-300
    % get indices of timepoints closest to -800 and -600 ms
    [~,il2] = min(abs(t + 800)); %-800
    [~,iu2] = min(abs(t + 600)); %-600

    % ==> choice
    cho = S.beh.choiceCat;
    % ==> context, contrast, orientation indicator variables
    ctx = S.exp.taskContext; ctr = S.exp.stimContrast; ori = S.exp.stimOriDeg;
    % x axis is orientation (the values differ for FN and JP!). Use unique values
    or = sort(unique(ori),'ascend')'; cx = sort(unique(ctx),'ascend')'; cr = sort(unique(ctr),'ascend')';            
    
    % ==> true dcCalAll (raw) values
    dvs_raw = S.dec.dvCatAll;
      
    for rc = 1:length(cr)
        for xc = 1:length(cx)
            % ==> index vector
            idv{iS,xc,rc} = (ori == 0) & (ctx == cx(xc)) & (ctr == cr(rc));
            
            % ==> store raw DV values
            dvsr{iS,xc,rc} = dvs_raw(idv{iS,xc,rc},:);
            
            % ==> store model fit DV values
            dvsm{iS,xc,rc} = dvs(idv{iS,xc,rc},:);
            
            % ==> store averages in two time windows
            dvw{iS,xc,rc} = [mean(mean(dvsm{iS,xc,rc}(:,il2:iu2),1),2), mean(mean(dvsm{iS,xc,rc}(:,il:iu),1),2)];
             
        end
    end
    
    mrks = {'-',':'}; clrs = {'r','b'};
    subplot(5,6,iS);
    hold on; hold all;    
    for rc = 1:length(cr)
        for xc = 1:length(cx)
            % ==> plot DV curves
            plot(t, mean(dvsm{iS,xc,rc}), [mrks{rc}, clrs{xc}]);
        end
    end
    ylim([-1,1]);
    drawnow; 
end

%%

% ==> F3A (J210906 low contrast, zero signal stimuli)
F3A();

%%
% ==> just DV graphs for low contrast (F3B)
F3B();

%%
% ==> make F3C histograms
F3C();

%% 
% ==> F3D scatterplots
F3D();

%%
% ==> eigenvector analysis
eig_test();