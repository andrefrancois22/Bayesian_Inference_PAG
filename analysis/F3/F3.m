clear all; close all; clc;
% ==> directories
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '../../data/';

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
% ==> context ccw (ctx == -1)
vccwlosF = [];    
% ==> context cw (ctx == 1)
vcwhisF = [];
% ==> context ccw (ctx == -1)
vccwhisF = [];     

% ==> store individual trial values
% ==> context cw (ctx == 1)
vcwlosJ = [];
% ==> context ccw (ctx == -1)
vccwlosJ = [];    
% ==> context cw (ctx == 1)
vcwhisJ = [];
% ==> context ccw (ctx == -1)
vccwhisJ = [];        

% % ==> store raw DV averages
% J_cw_lo_raw_DVs = [];
% J_cw_hi_raw_DVs = [];
% J_ccw_lo_raw_DVs = [];
% J_ccw_hi_raw_DVs = [];
% 
% % ==> store raw DV averages
% F_cw_lo_raw_DVs = [];
% F_cw_hi_raw_DVs = [];
% F_ccw_lo_raw_DVs = [];
% F_ccw_hi_raw_DVs = [];

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
        
%         % ==> store averages of raw DVs 
%         F_cw_lo_raw_DVs(iSfln,:)  = mean(dvs_raw(idxcwlo,:),1);
%         F_cw_hi_raw_DVs(iSfln,:)  = mean(dvs_raw(idxcwhi,:),1);
%         F_ccw_lo_raw_DVs(iSfln,:) = mean(dvs_raw(idxccwlo,:),1);
%         F_ccw_hi_raw_DVs(iSfln,:) = mean(dvs_raw(idxccwhi,:),1);        
        
        % ==> store all single trial dv (model fit) initial offsets 
        % normalized by average dynamic range
        vcwlosF  = [vcwlosF;  vcwlo];
        vccwlosF = [vccwlosF; vccwlo];
        vcwhisF  = [vcwhisF;  vcwhi];
        vccwhisF = [vccwhisF; vccwhi];      
        
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
        
%         % ==> store averages of raw DVs 
%         J_cw_lo_raw_DVs(iSfln  -13,:)  = mean(dvs_raw(idxcwlo,:),1);
%         J_cw_hi_raw_DVs(iSfln  -13,:)  = mean(dvs_raw(idxcwhi,:),1);
%         J_ccw_lo_raw_DVs(iSfln  -13,:) = mean(dvs_raw(idxccwlo,:),1);
%         J_ccw_hi_raw_DVs(iSfln  -13,:) = mean(dvs_raw(idxccwhi,:),1);
        
        % ==> store all single trial dv (model fit) initial offsets 
        % normalized by average dynamic range
        vcwlosJ  = [vcwlosJ;  vcwlo];
        vccwlosJ = [vccwlosJ; vccwlo];
        vcwhisJ  = [vcwhisJ;  vcwhi];
        vccwhisJ = [vccwhisJ; vccwhi];     
        
    end    
end

%%
% ==> histograms for F3

% ==> make F3C histograms
F3C();

%%
% ==> just DV graphs for low contrast (F3B)
F3B();

%%
% ==> F3D scatterplots
F3D();

%% 
% ==> F3A (J210906 low contrast, zero signal stimuli)
F3A();
    
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
% ==> bootstrap sample number
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