clear all; close all; clc;
% ==> directories
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
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
% ==> histograms for F3

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
    
%% ==> estimate the slope of the line in F3 scatterplot
       
F = vertcat(dvw{1:13,:,:});
J = vertcat(dvw{14:29,:,:});

FJ = [J;F];

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
figure(); clf; set(gcf,'Color','w');
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