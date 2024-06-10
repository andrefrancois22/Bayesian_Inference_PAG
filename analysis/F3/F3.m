clear all; close all; clc;
% ==> directories
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '../data/';

%%
close all; clc;

figure(); set(gcf,'Color','w');

% => cw cases
Fmu_cw_lo = nan(13,2);
Fmu_cw_hi = nan(13,2);

Jmu_cw_lo = nan(16,2);
Jmu_cw_hi = nan(16,2);

% => ccw cases
Fmu_ccw_lo = nan(13,2);
Fmu_ccw_hi = nan(13,2);

Jmu_ccw_lo = nan(16,2);
Jmu_ccw_hi = nan(16,2);

% ==> for storing all neutral DV trajectories
% => cw cases
F_cw_lo = nan(13,200);
F_cw_hi = nan(13,200);
J_cw_lo = nan(16,200);
J_cw_hi = nan(16,200);

% => ccw cases
F_ccw_lo = nan(13,200);
F_ccw_hi = nan(13,200);
J_ccw_lo = nan(16,200);
J_ccw_hi = nan(16,200);

% ==> store individual trial values
% ==> context cw (ctx == 1)
vcwlosF = [];
% ==> context ccw (ctx == 1)
vccwlosF = [];    
% ==> context cw (ctx == 1)
vcwhisF = [];
% ==> context ccw (ctx == 1)
vccwhisF = [];     

% ==> store individual trial values
% ==> context cw (ctx == 1)
vcwlosJ = [];
% ==> context ccw (ctx == 1)
vccwlosJ = [];    
% ==> context cw (ctx == 1)
vcwhisJ = [];
% ==> context ccw (ctx == 1)
vccwhisJ = [];        

% ==> store raw DV averages
J_cw_lo_raw_DVs = [];
J_cw_hi_raw_DVs = [];
J_ccw_lo_raw_DVs = [];
J_ccw_hi_raw_DVs = [];

% ==> store raw DV averages
F_cw_lo_raw_DVs = [];
F_cw_hi_raw_DVs = [];
F_ccw_lo_raw_DVs = [];
F_ccw_hi_raw_DVs = [];

% ==> store dynamic range means
mus_dynr = nan(1,29);

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
    or = unique(ori)'; cx = unique(ctx)'; cr = unique(ctr)';
    
    % ==> dynamic range for session
    dynr = max([max(dvs, [], 2) - dvs(:,1), -(min(dvs, [], 2) - dvs(:,1))], [], 2);
    % ==> average of the dynamic range for that session
    dynr_mu = mean(dynr);
    % ==> store the average for that session
    mus_dynr(iSfln) = dynr_mu;
    
    % indices
    idxcwlo =  (ori == 0) & (ctx == cx(2)) & (ctr == cr(1));
    idxccwlo = (ori == 0) & (ctx == cx(1)) & (ctr == cr(1));
    
    idxcwhi =  (ori == 0) & (ctx == cx(2)) & (ctr == cr(2));
    idxccwhi = (ori == 0) & (ctx == cx(1)) & (ctr == cr(2));    
    
    % ==> context cw (ctx == 1)
    vcwlo  = dvs(idxcwlo,:);
    % ==> context ccw (ctx == 1)
    vccwlo = dvs(idxccwlo,:);    
    % ==> context cw (ctx == 1)
    vcwhi  = dvs(idxcwhi,:);
    % ==> context ccw (ctx == 1)
    vccwhi = dvs(idxccwhi,:);                   
    
    % ==> true dcCalAll (raw) values
    dvs_raw = S.dec.dvCatAll;
      
    subplot(5,6,iSfln);
    hold on; hold all;    
    plot(t,mean(vcwlo,1),  'r:');
    plot(t,mean(vccwlo,1), 'b:');
    plot(t,mean(vcwhi,1),  'r-');
    plot(t,mean(vccwhi,1), 'b-');    
    ylim([-1,1]);
    drawnow;

    % => if Monkey F
    if iSfln <= 13
        % ==> store averages in two time windows
        % => lo contrast       
        Fmu_cw_lo(iSfln,:)  = [mean(mean(vcwlo(:,il2:iu2),1),2), mean(mean(vcwlo(:,il:iu),1),2)];
        % => hi contrast
        Fmu_cw_hi(iSfln,:)  = [mean(mean(vcwhi(:,il2:iu2),1),2), mean(mean(vcwhi(:,il:iu),1),2)];        
        % => lo contrast
        Fmu_ccw_lo(iSfln,:) = [mean(mean(vccwlo(:,il2:iu2),1),2), mean(mean(vccwlo(:,il:iu),1),2)];
        % => hi contrast
        Fmu_ccw_hi(iSfln,:) = [mean(mean(vccwhi(:,il2:iu2),1),2), mean(mean(vccwhi(:,il:iu),1),2)];    
        
        % ==> store average DV trajectories        
        F_cw_lo(iSfln,:)  = mean(vcwlo,1);    
        F_cw_hi(iSfln,:)  = mean(vcwhi,1);  
        F_ccw_lo(iSfln,:) = mean(vccwlo,1);    
        F_ccw_hi(iSfln,:) = mean(vccwhi,1);  
        
        % ==> store averages of raw DVs 
        F_cw_lo_raw_DVs(iSfln,:)  = mean(dvs_raw(idxcwlo,:),1);
        F_cw_hi_raw_DVs(iSfln,:)  = mean(dvs_raw(idxcwhi,:),1);
        F_ccw_lo_raw_DVs(iSfln,:) = mean(dvs_raw(idxccwlo,:),1);
        F_ccw_hi_raw_DVs(iSfln,:) = mean(dvs_raw(idxccwhi,:),1);        
        
        % ==> store all single trial dv (model fit) initial offsets 
        % normalized by average dynamic range
        vcwlosF  = [vcwlosF;  vcwlo];%  / dynr_mu ];
        vccwlosF = [vccwlosF; vccwlo];% / dynr_mu ]; 
        vcwhisF  = [vcwhisF;  vcwhi];%  / dynr_mu ];
        vccwhisF = [vccwhisF; vccwhi];% / dynr_mu ];      
        
    % => if Monkey J    
    elseif iSfln >= 14
        % ==> store averages in two time windows        
        % => lo contrast
        Jmu_cw_lo(iSfln -13,:)  = [mean(mean(vcwlo(:,il2:iu2),1),2), mean(mean(vcwlo(:,il:iu),1),2)];
        % => hi contrast
        Jmu_cw_hi(iSfln -13,:)  = [mean(mean(vcwhi(:,il2:iu2),1),2), mean(mean(vcwhi(:,il:iu),1),2)];        
        % => lo contrast
        Jmu_ccw_lo(iSfln -13,:) = [mean(mean(vccwlo(:,il2:iu2),1),2), mean(mean(vccwlo(:,il:iu),1),2)];
        % => hi contrast
        Jmu_ccw_hi(iSfln -13,:) = [mean(mean(vccwhi(:,il2:iu2),1),2), mean(mean(vccwhi(:,il:iu),1),2)];      
        
        % ==> store average DV trajectories (model fits)       
        J_cw_lo(iSfln  -13,:) = mean(vcwlo,1);    
        J_cw_hi(iSfln  -13,:) = mean(vcwhi,1);  
        J_ccw_lo(iSfln -13,:) = mean(vccwlo,1);    
        J_ccw_hi(iSfln -13,:) = mean(vccwhi,1);    
        
        % ==> store averages of raw DVs 
        J_cw_lo_raw_DVs(iSfln  -13,:)  = mean(dvs_raw(idxcwlo,:),1);
        J_cw_hi_raw_DVs(iSfln  -13,:)  = mean(dvs_raw(idxcwhi,:),1);
        J_ccw_lo_raw_DVs(iSfln  -13,:) = mean(dvs_raw(idxccwlo,:),1);
        J_ccw_hi_raw_DVs(iSfln  -13,:) = mean(dvs_raw(idxccwhi,:),1);
        
        % ==> store all single trial dv (model fit) initial offsets 
        % normalized by average dynamic range
        vcwlosJ  = [vcwlosJ;  vcwlo];%  / dynr_mu ];
        vccwlosJ = [vccwlosJ; vccwlo];% / dynr_mu ]; 
        vcwhisJ  = [vcwhisJ;  vcwhi];%  / dynr_mu ];
        vccwhisJ = [vccwhisJ; vccwhi];% / dynr_mu ];      
        
    end    
end

%% ==> histograms for F3
close all; clc;

% ==> all CW offset values for vertical stimulus - across hi and lo contrasts
cwsF  = [vcwlosF; vcwhisF];   
% ==> initial offset
cwsFi  = cwsF(:,1);
% ==> all CCW offset values for vertical stimulus - across hi and lo contrasts
ccwsF = [vccwlosF; vccwhisF]; 
% ==> initial offset
ccwsFi = ccwsF(:,1);

% ==> all CW offset values for vertical stimulus - across hi and lo contrasts
cwsJ  = [vcwlosJ; vcwhisJ];   
% ==> initial offset
cwsJi  = cwsJ(:,1);
% ==> all CCW offset values for vertical stimulus - across hi and lo contrasts
ccwsJ = [vccwlosJ; vccwhisJ]; 
% ==> initial offset
ccwsJi = ccwsJ(:,1);

% ==>
cwsF_mu  = mean(cwsF(:,il2:iu2),2);
ccwsF_mu = mean(ccwsF(:,il2:iu2),2);

cwsJ_mu  = mean(cwsJ(:,il2:iu2),2);
ccwsJ_mu = mean(ccwsJ(:,il2:iu2),2);

% ==> edges
edgs = -1.5:0.1:1.5;

fg = figure(); set(fg,'color','white'); set(fg, 'Position',[675 421 660 541])
subplot(2,2,1);
hFcw = histogram(cwsF_mu,edgs,'facecolor',[0.15,0.75,0.5]); hold on; hold all;
plot(zeros(1,max(hFcw.Values)),1:max(hFcw.Values),'r--','linewidth',2)
title('Animal F: cw prior context')
xlabel('Signed DV initial offset')
ylabel('Frequency')
drawnow;
subplot(2,2,2);
hFccw = histogram(ccwsF_mu,edgs,'facecolor',[0.15,0.75,0.5]); hold on; hold all;
plot(zeros(1,max(hFccw.Values)),1:max(hFccw.Values),'r--','linewidth',2)
title('Animal F: ccw prior context')
xlabel('Signed DV initial offset')
ylabel('Frequency')
drawnow;
 
subplot(2,2,3);
hJcw = histogram(cwsJ_mu,edgs,'facecolor',[1,0.5,0]); hold on; hold all;
plot(zeros(1,max(hJcw.Values)),1:max(hJcw.Values),'r--','linewidth',2)
title('Animal J: cw prior context')
xlabel('Signed DV initial offset')
ylabel('Frequency')
drawnow;
subplot(2,2,4);
hJccw = histogram(ccwsJ_mu,edgs,'facecolor',[1,0.5,0]); hold on; hold all;
plot(zeros(1,max(hJccw.Values)),1:max(hJccw.Values),'r--','linewidth',2)
title('Animal J: ccw prior context')
xlabel('Signed DV initial offset')
ylabel('Frequency')
drawnow;

% ==> 
% ==> signrank wilcoxon test
[p_cwJ,~]  = signrank(cwsJ_mu);
[p_ccwJ,~] = signrank(ccwsJ_mu);

[p_cwF,~]  = signrank(cwsF_mu);
[p_ccwF,~] = signrank(ccwsF_mu);

fprintf('Wilcoxon test p values evaluating initial offsets in the vertical stimulus condition...\n')
fprintf('Wilcoxon test for animal F cw cases, p = %d...\n', p_cwF)
fprintf('Wilcoxon test for animal F ccw cases, p = %d...\n', p_ccwF)


fprintf('Wilcoxon test for animal J cw cases, p = %d...\n', p_cwJ)
fprintf('Wilcoxon test for animal J ccw cases, p = %d...\n', p_ccwJ)
 
% ==> ranksum
[PJ,~] = ranksum(cwsJ_mu,ccwsJ_mu)

[PF,~] = ranksum(cwsF_mu,ccwsF_mu)

% ==> restrict to the low contrast trials and repeat ranksum


%%
% close all; clc;

figure(); set(gcf,'Color','w'); set(gcf,'Position',[86 490 1834 472]);
subplot(1,3,1);
hold on; hold all;
% => monkey F
s1 = scatter(Fmu_cw_lo(:,1), Fmu_cw_lo(:,2),  90, 'o','markerfacecolor', [0.25,0.85,0.6], 'markeredgecolor', 'r'); %, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s2 = scatter(Fmu_ccw_lo(:,1),Fmu_ccw_lo(:,2), 90, 'o','markerfacecolor', [0.25,0.85,0.6], 'markeredgecolor', 'b'); %, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s3 = scatter(Fmu_cw_hi(:,1), Fmu_cw_hi(:,2),  90, 'o','markerfacecolor', [0.15,0.75,0.5], 'markeredgecolor', 'r');
s4 = scatter(Fmu_ccw_hi(:,1),Fmu_ccw_hi(:,2), 90, 'o','markerfacecolor', [0.15,0.75,0.5], 'markeredgecolor', 'b');
% => monkey J
s5 = scatter(Jmu_cw_lo(:,1), Jmu_cw_lo(:,2),  90, 'o','markerfacecolor', [1,0.75,0], 'markeredgecolor', 'r'); %, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s6 = scatter(Jmu_ccw_lo(:,1),Jmu_ccw_lo(:,2), 90, 'o','markerfacecolor', [1,0.75,0], 'markeredgecolor', 'b'); %, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s7 = scatter(Jmu_cw_hi(:,1), Jmu_cw_hi(:,2),  90, 'o','markerfacecolor', [1,0.5,0], 'markeredgecolor', 'r');
s8 = scatter(Jmu_ccw_hi(:,1),Jmu_ccw_hi(:,2), 90, 'o','markerfacecolor', [1,0.5,0], 'markeredgecolor', 'b');
hold on; hold all
plot(-0.6:0.1:0.6,-0.6:0.1:0.6,'k--')
plot(-0.6:0.1:0.6,zeros([1,length(-0.6:0.1:0.6)]),'k--')
plot(zeros([1,length(-0.6:0.1:0.6)]),-0.6:0.1:0.6,'k--')
axis equal; axis square;
xlabel('Average (-800ms to -600ms) window')
ylabel('Average (-500ms to -300ms) window')
legend([s1,s2,s3,s4,s5,s6,s7,s8], 'Monkey F cw low contrast','Monkey F ccw low contrast','Monkey F cw high contrast','Monkey F ccw high contrast', ...
                                  'Monkey J cw low contrast','Monkey J ccw low contrast','Monkey J cw high contrast','Monkey J ccw high contrast',... 
                                 'Location','SouthEast')

subplot(1,3,2);
hold on; hold all;
p1 = plot(linspace(-790,-40,200), mean(J_cw_lo,1), 'color', [1,0,0,0.35],  'linewidth',8);
p2 = plot(linspace(-790,-40,200), mean(J_cw_hi,1), 'r-',  'linewidth',4);
p3 = plot(linspace(-790,-40,200), mean(J_ccw_lo,1), 'color', [0,0,1,0.35], 'linewidth',8);
p4 = plot(linspace(-790,-40,200), mean(J_ccw_hi,1), 'b-', 'linewidth',4);
% => time windows
plot(repmat(-800,[20,1]),linspace(-0.2,0.2,20),'k--')
plot(repmat(-600,[20,1]),linspace(-0.2,0.2,20),'k--')

plot(repmat(-500,[20,1]),linspace(-0.2,0.2,20),'k:')
plot(repmat(-300,[20,1]),linspace(-0.2,0.2,20),'k:')
xlim([-840,-30])
title('Monkey J');
xlabel('Time');
ylabel('Signed DV average (vertical stimulus)');
legend([p1,p2,p3,p4],'cw context, low contrast', 'cw context high contrast', ...
                     'ccw context, low contrast', 'ccw context, high contrast');

subplot(1,3,3);
hold on; hold all;
p1 = plot(linspace(-795,-45,200), mean(F_cw_lo,1), 'color', [1,0,0,0.35],  'linewidth',8);
p2 = plot(linspace(-795,-45,200), mean(F_cw_hi,1), 'r-',  'linewidth',4);
p3 = plot(linspace(-795,-45,200), mean(F_ccw_lo,1), 'color', [0,0,1,0.35], 'linewidth',8);
p4 = plot(linspace(-795,-45,200), mean(F_ccw_hi,1), 'b-', 'linewidth',4);
% => time windows
plot(repmat(-800,[20,1]),linspace(-0.2,0.2,20),'k--')
plot(repmat(-600,[20,1]),linspace(-0.2,0.2,20),'k--')

plot(repmat(-500,[20,1]),linspace(-0.2,0.2,20),'k:')
plot(repmat(-300,[20,1]),linspace(-0.2,0.2,20),'k:')
xlim([-840,-30])
title('Monkey F');
xlabel('Time');
ylabel('Signed DV average (vertical stimulus)');
legend([p1,p2,p3,p4],'cw context, low contrast', 'cw context high contrast', ...
                     'ccw context, low contrast', 'ccw context, high contrast');


%% ==> just DV graphs for low contrast

figure();  set(gcf,'Color','w'); set(gcf, 'Position',[675 553 963 409])
subplot(1,2,1);
hold on; hold all;
p1 = plot(linspace(-790,-40,200), mean(J_cw_lo,1), 'color', [1,0,0,0.35],  'linewidth',8);
p3 = plot(linspace(-790,-40,200), mean(J_ccw_lo,1), 'color', [0,0,1,0.35], 'linewidth',8);
% ==> timeSac(16) == -790 and timeSac(31) == -40
p4 = plot(timeSac(16:31), mean(J_cw_lo_raw_DVs(:,16:31),1), 'o', 'color', [1,0,0,0.35]);
p5 = plot(timeSac(16:31), mean(J_ccw_lo_raw_DVs(:,16:31),1), 'o', 'color', [0,0,1,0.35], 'linewidth',1);
% => time windows
plot(repmat(-800,[20,1]),linspace(-0.2,0.2,20),'k--')
plot(repmat(-600,[20,1]),linspace(-0.2,0.2,20),'k--')

plot(repmat(-500,[20,1]),linspace(-0.2,0.2,20),'k:')
plot(repmat(-300,[20,1]),linspace(-0.2,0.2,20),'k:')
xlim([-840,-30])
title('Monkey J');
xlabel('Time');
ylabel('Signed DV average (vertical stimulus)');
% legend([p1,p3],'cw context, low contrast', ...
%                'ccw context, low contrast');

subplot(1,2,2);
hold on; hold all;
p1 = plot(linspace(-795,-45,200), mean(F_cw_lo,1), 'color', [1,0,0,0.35],  'linewidth',8);
p3 = plot(linspace(-795,-45,200), mean(F_ccw_lo,1), 'color', [0,0,1,0.35], 'linewidth',8);
% ==> timeSac(16) == -790 and timeSac(31) == -40
p4 = plot(timeSac(16:31), mean(F_cw_lo_raw_DVs(:,16:31),1), 'o', 'color', [1,0,0,0.35]);
p5 = plot(timeSac(16:31), mean(F_ccw_lo_raw_DVs(:,16:31),1), 'o', 'color', [0,0,1,0.35], 'linewidth',1);
% => time windows
plot(repmat(-800,[20,1]),linspace(-0.2,0.2,20),'k--')
plot(repmat(-600,[20,1]),linspace(-0.2,0.2,20),'k--')

plot(repmat(-500,[20,1]),linspace(-0.2,0.2,20),'k:')
plot(repmat(-300,[20,1]),linspace(-0.2,0.2,20),'k:')
xlim([-840,-30])
title('Monkey F');
xlabel('Time');
ylabel('Signed DV average (vertical stimulus)');
legend([p1,p3],'cw context, low contrast', ...
               'ccw context, low contrast');
                 
%% ==> just seprate out by monkey

figure(); set(gcf,'Color','w'); set(gcf,'Position',[178 358 1325 604]);
subplot(1,2,1);
hold on; hold all;
% => monkey F
s1 = scatter(Fmu_cw_lo(:,1), Fmu_cw_lo(:,2),  140, 'o','markerfacecolor', [0.25,0.85,0.6], 'markeredgecolor', 'r');%, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s2 = scatter(Fmu_ccw_lo(:,1),Fmu_ccw_lo(:,2), 140, 'o','markerfacecolor', [0.25,0.85,0.6], 'markeredgecolor', 'b');%, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s3 = scatter(Fmu_cw_hi(:,1), Fmu_cw_hi(:,2),  140, 'o','markerfacecolor', [0.15,0.75,0.5], 'markeredgecolor', 'r');
s4 = scatter(Fmu_ccw_hi(:,1),Fmu_ccw_hi(:,2), 140, 'o','markerfacecolor', [0.15,0.75,0.5], 'markeredgecolor', 'b');
% ==> lines
plot(-0.6:0.1:0.6,-0.6:0.1:0.6,'k--')
plot(-0.6:0.1:0.6,zeros([1,length(-0.6:0.1:0.6)]),'k--')
plot(zeros([1,length(-0.6:0.1:0.6)]),-0.6:0.1:0.6,'k--')
axis equal; axis square;
xlabel('Average (-800ms to -600ms) window');
ylabel('Average (-500ms to -300ms) window');
legend([s1,s2,s3,s4],'cw low contrast','ccw low contrast','cw high contrast','ccw high contrast','Location','SouthEast')
title('Monkey F');

subplot(1,2,2);
hold on; hold all;
% => monkey J
s1 = scatter(Jmu_cw_lo(:,1), Jmu_cw_lo(:,2),  140, 'o','markerfacecolor', [1,0.75,0], 'markeredgecolor', 'r');%, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s2 = scatter(Jmu_ccw_lo(:,1),Jmu_ccw_lo(:,2), 140, 'o','markerfacecolor', [1,0.75,0], 'markeredgecolor', 'b');%, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s3 = scatter(Jmu_cw_hi(:,1), Jmu_cw_hi(:,2),  140, 'o','markerfacecolor', [1,0.5,0], 'markeredgecolor', 'r');
s4 = scatter(Jmu_ccw_hi(:,1),Jmu_ccw_hi(:,2), 140, 'o','markerfacecolor', [1,0.5,0], 'markeredgecolor', 'b');
% ==> lines
plot(-0.6:0.1:0.6,-0.6:0.1:0.6,'k--')
plot(-0.6:0.1:0.6,zeros([1,length(-0.6:0.1:0.6)]),'k--')
plot(zeros([1,length(-0.6:0.1:0.6)]),-0.6:0.1:0.6,'k--')
axis equal; axis square;
xlabel('Average (-800ms to -600ms) window')
ylabel('Average (-500ms to -300ms) window')
legend([s1,s2,s3,s4],'cw low contrast','ccw low contrast','cw high contrast','ccw high contrast','Location','SouthEast')
title('Monkey J');

% text(Jmu_cw_lo(:,1)  - 0.015, Jmu_cw_lo(:,2),  string([14:29]),'fontsize',6)
% text(Jmu_ccw_lo(:,1) - 0.015,Jmu_ccw_lo(:,2), string([14:29]),'fontsize',6)
% text(Jmu_cw_hi(:,1)  - 0.015, Jmu_cw_hi(:,2),  string([14:29]),'fontsize',6)
% text(Jmu_ccw_hi(:,1) - 0.015,Jmu_ccw_hi(:,2), string([14:29]),'fontsize',6)


%% ==> just one example
% close all; clc;

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

% hold on; hold all;   
% p3 = plot(t,mean(vcwhi,1),  'r-','linewidth',4);
% p4 = plot(t,mean(vccwhi,1), 'b-','linewidth',4);    

% legend([p1,p2,p3,p4],'CW prior low contrast','CCW prior low contrast', 'CW prior high contrast','CCW prior high contrast');
legend([p1,p2],'cw context, low contrast', ...
               'ccw context, low contrast');
drawnow;
    
%% ==> estimate the slope of the line in F3 scatterplot

% ==> Noisy x and y variables - use first eigenvector

F = [Fmu_cw_lo; Fmu_cw_hi; Fmu_ccw_lo; Fmu_ccw_hi];
J = [Jmu_cw_lo; Jmu_cw_hi; Jmu_ccw_lo; Jmu_ccw_hi];
FJ = [F;J];

% ==> covariance matrices
covF = cov(F);
covJ = cov(J);
% ==> covariance matrix for all the data
covFJ = cov(FJ);

% ==> eigendecomposition
[V,D] = eig(covFJ);

% ==> get index of max eigenvalue
[~,I] = max(D);

% ==> mean of the distribution
muFJ = mean(FJ)';

% ==> angle of first eigenvector
[TH,~] = cart2pol(V(1,I(2)),V(2,I(2)));

% ==> visualize first eigenvector over distribution
figure(1); clf; set(gcf,'Color','w');
scatter(FJ(:,1),FJ(:,2),'ko','filled');
hold on; hold all;
plot([muFJ(1),V(1,I(2))],[muFJ(2),V(2,I(2))],'r-')
% ==> lines
plot(-0.6:0.1:0.6,-0.6:0.1:0.6,'k--')
plot(-0.6:0.1:0.6,zeros([1,length(-0.6:0.1:0.6)]),'k--')
plot(zeros([1,length(-0.6:0.1:0.6)]),-0.6:0.1:0.6,'k--')
axis equal; axis square;
xlabel('Average (-800ms to -600ms) window');
ylabel('Average (-500ms to -300ms) window');
title(['First PC slope (rad): ', num2str(TH)])

%% ==> random sampling

% ==> Noisy x and y variables - use first eigenvector

F = [Fmu_cw_lo; Fmu_cw_hi; Fmu_ccw_lo; Fmu_ccw_hi];
J = [Jmu_cw_lo; Jmu_cw_hi; Jmu_ccw_lo; Jmu_ccw_hi];
FJ = [F;J];

% ==> covariance matrices
covF = cov(F);
covJ = cov(J);
% ==> covariance matrix for all the data

% ==> for random sample with replacement
N = length(FJ);
K = length(FJ);
% ==> 100 bootstrap samples
BT = 1000;

% ==> store slopes
s = nan(1,BT);

% ==> sample
for b = 1:BT

    % indices for random sample with replacement
    idx = randsample(N,K,true);
    FJs = FJ(idx,:);

    covFJs = cov(FJs);

    % ==> eigendecomposition
    [V,D] = eig(covFJs);

    % ==> get index of max eigenvalue
    [~,I] = max(D);

    % ==> mean of the distribution
    muFJ = mean(FJs)';

    % ==> angle of first eigenvector
    [TH,~] = cart2pol(V(1,I(2)),V(2,I(2)));

    % ==> store slope
    s(b) = TH;
    
    % ==> visualize first eigenvector over distribution
    fg = figure(1); fg.Position = [333 552 850 338]; 
    clf; set(gcf,'Color','w');
    subplot(1,2,1);
    scatter(FJs(:,1),FJs(:,2),'ko','filled');
    hold on; hold all;
    plot([muFJ(1),V(1,I(2))],[muFJ(2),V(2,I(2))],'r-'); 
    % ==> lines
    plot(-0.6:0.1:0.6,-0.6:0.1:0.6,'r--')
    plot(-0.6:0.1:0.6,zeros([1,length(-0.6:0.1:0.6)]),'k--')
    plot(zeros([1,length(-0.6:0.1:0.6)]),-0.6:0.1:0.6,'k--')
    axis equal; axis square;
    xlabel('Average (-800ms to -600ms) window');
    ylabel('Average (-500ms to -300ms) window');
    title(['First PC slope (rad): ', num2str(TH)])
    drawnow;

    % ==> update in terminal
    fprintf('computed first PC slope for sampled data %d of %d...\n',b,BT)
end

% ==> eigendecomposition
[V,D] = eig(covFJs);
% ==> get index of max eigenvalue
[~,I] = max(D);
% ==> mean of the distribution
muFJ = mean(FJ)';

% ==> angle of first eigenvector
[TH,~] = cart2pol(V(1,I(2)),V(2,I(2)));
% ==> visualize first eigenvector over distribution
fg = figure(1); fg.Position = [333 552 850 338]; 
clf; set(gcf,'Color','w');
subplot(1,2,1);
scatter(FJ(:,1),FJ(:,2),'ko','filled');
hold on; hold all;
plot([muFJ(1),V(1,I(2))],[muFJ(2),V(2,I(2))],'r-'); 
% ==> lines
plot(-0.6:0.1:0.6,-0.6:0.1:0.6,'r--')
plot(-0.6:0.1:0.6,zeros([1,length(-0.6:0.1:0.6)]),'k--')
plot(zeros([1,length(-0.6:0.1:0.6)]),-0.6:0.1:0.6,'k--')
axis equal; axis square;
xlabel('Average (-800ms to -600ms) window');
ylabel('Average (-500ms to -300ms) window');
title(['First PC slope (rad): ', num2str(TH)])
drawnow;

bns = 25;
subplot(1,2,2); cla;
h = histogram(s,bns,'Facecolor','r'); hold on; hold all;
xlabel('First PC slope (radians)');
ylabel('frequency');
title('distribution of slopes of first PC (JP)');
plot(repmat(pi/4,[1,length(0:max(h.Values))]),0:max(h.Values),'r--')
text(pi/4 + 0.025,10,'\pi/4 equality line (x=y) angle (slope)','Rotation',90,'Color','r')
xlim([pi/5,pi/2]);

% ==> compute difference between estimated slopes and 'null' slope (e.g. x=y, or TH = pi/4)
ds = s - pi/4;

% ==> test that ds deviated significantly from zero (means slope angle is greater than pi/4 e.g. x=y)
[ps, ~] = signrank(ds)