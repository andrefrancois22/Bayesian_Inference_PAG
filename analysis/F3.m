clear all; close all; clc;
% ==> directories
dataPath     = strcat('/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data');
drc = '../data/';

%% ==> accumulation bound panel

figure(3); set(gcf,'Color','w'); set(gcf, 'Position',[675 461 932 501])

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for fc = 1:0.1:5

pr_cw  =  0.5 * fc;
pr_ccw = -0.5 * fc;

% => number of trials
N = 500;
% => number of timepoints
tm = 200;

% ==>
sens  = randn(N,tm);
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

% figure(1); set(gcf,'Color','w'); set(gcf,'Position',[124 717 1682 245]);
% subplot(1,5,1); 
% hold on; hold all;
% plot((csens + pr_cw)', 'color', [1,0.5,0.5,0.5]);     ylim([-bnd - 10, bnd + 10]);
% plot(1:t,repmat( bnd,tm), 'r--', 'linewidth', 4);
% plot(1:t,repmat(-bnd,tm), 'b--', 'linewidth', 4);
% title('cw prior ctx')
% subplot(1,5,2); 
% hold on; hold all;
% plot((csens + pr_ccw)', 'color', [0.5,0.5,1,0.5]);     ylim([-bnd - 10, bnd + 10]);
% plot(1:t,repmat( bnd,tm), 'r--', 'linewidth', 4);
% plot(1:t,repmat(-bnd,tm), 'b--', 'linewidth', 4);
% title('ccw prior ctx')
%  
% subplot(1,5,3); 
% hold on; hold all;
% plot(1:t,repmat( bnd,tm), 'r--', 'linewidth', 4);
% plot(1:t,repmat(-bnd,tm), 'b--', 'linewidth', 4);
% plot(csensb_cw', 'color', [1,0.5,0.5,0.5]);     ylim([-bnd - 10, bnd + 10]);
% title('cw prior ctx')
% subplot(1,5,4); 
% hold on; hold all;
% plot(1:t,repmat( bnd,tm), 'r--', 'linewidth', 4);
% plot(1:t,repmat(-bnd,tm), 'b--', 'linewidth', 4);
% plot(csensb_ccw', 'color', [0.5,0.5,1,0.5]);     ylim([-bnd - 10, bnd + 10]);
% title('ccw prior ctx')
% 
% subplot(1,5,5);
% plot(mean(csensb_cw,1),'linewidth',1, 'color',  [1,0.5,0.5,0.5]);
% hold on; hold all;
% plot(mean(csensb_ccw,1),'linewidth',1, 'color', [0.5,0.5,1,0.5]);
% ylabel('Average signed DV (acc. bound)')
% xlabel('Time')
% title('Average over trials')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%==> case with a bias that grows linearly with t (mean drift)

pr_cw_d =   fc*(0:(tm-1))/(0.5*tm-1);
pr_ccw_d = -fc*(0:(tm-1))/(0.5*tm-1);

% ==> check that average bias is equal to pr_cw and pr_ccw
% assert(mean(pr_cw_d)  == pr_cw);
% assert(mean(pr_ccw_d) == pr_ccw);

% ==> csens with bound - two cases: with ccw prior or with cw prior
csensb_cw_d =  csens + repmat(pr_cw_d, [N,1]);
csensb_ccw_d = csens + repmat(pr_ccw_d,[N,1]);

% ==> set DV values for timepoints over trials that reach or exceed a bound
% to the bound value.
% ==> cw prior ctx
csensb_cw_d(csensb_cw_d >=  bnd) =  bnd;
csensb_cw_d(csensb_cw_d <= -bnd) = -bnd;
% ==> ccw prior ctx
csensb_ccw_d(csensb_ccw_d >=  bnd) =  bnd;
csensb_ccw_d(csensb_ccw_d <= -bnd) = -bnd;

% ==> once a bnd has been reached at time t for a given trial, set the remaining DV values
% starting from t+1 to the bound for that trial
% ==> cw cases
for i = 1:N
    bidx = min(find(csensb_cw_d(i,:) >=  bnd));
    % => set to bound
    if ~isempty(bidx)
        csensb_cw_d(i,bidx:end) = bnd;
    end
end
for i = 1:N
    bidx = min(find(csensb_cw_d(i,:) <=  -bnd));
    % => set to bound
    if ~isempty(bidx)
        csensb_cw_d(i,bidx:end) = -bnd;
    end
end
% ==> ccw context cases
for i = 1:N
    bidx = min(find(csensb_ccw_d(i,:) >=  bnd));
    % => set to bound
    if ~isempty(bidx)
        csensb_ccw_d(i,bidx:end) = bnd;
    end
end
for i = 1:N
    bidx = min(find(csensb_ccw_d(i,:) <=  -bnd));
    % => set to bound
    if ~isempty(bidx)
        csensb_ccw_d(i,bidx:end) = -bnd;
    end
end

% figure(2); set(gcf,'Color','w'); set(gcf,'Position',[124 717 1682 245]);
% subplot(1,5,1); 
% hold on; hold all;
% plot((csens + repmat(pr_cw_d, [N,1]))', 'color', [1,0.5,0.5,0.5]);     ylim([-bnd - 10, bnd + 10]);
% plot(1:t,repmat( bnd,tm), 'r--', 'linewidth', 4);
% plot(1:t,repmat(-bnd,tm), 'b--', 'linewidth', 4);
% title('cw prior ctx')
% subplot(1,5,2); 
% hold on; hold all;
% plot((csens + repmat(pr_ccw_d, [N,1]))', 'color', [0.5,0.5,1,0.5]);     ylim([-bnd - 10, bnd + 10]);
% plot(1:t,repmat( bnd,tm), 'r--', 'linewidth', 4);
% plot(1:t,repmat(-bnd,tm), 'b--', 'linewidth', 4);
% title('ccw prior ctx')
%  
% subplot(1,5,3); 
% hold on; hold all;
% plot(1:t,repmat( bnd,tm), 'r--', 'linewidth', 4);
% plot(1:t,repmat(-bnd,tm), 'b--', 'linewidth', 4);
% plot(csensb_cw_d', 'color', [1,0.5,0.5,0.5]);     ylim([-bnd - 10, bnd + 10]);
% title('cw prior ctx')
% subplot(1,5,4); 
% hold on; hold all;
% plot(1:t,repmat( bnd,tm), 'r--', 'linewidth', 4);
% plot(1:t,repmat(-bnd,tm), 'b--', 'linewidth', 4);
% plot(csensb_ccw_d', 'color', [0.5,0.5,1,0.5]);     ylim([-bnd - 10, bnd + 10]);
% title('ccw prior ctx')
% 
% subplot(1,5,5);
% plot(mean(csensb_cw_d,1),'linewidth',1, 'color',  [1,0.5,0.5,0.5]);
% hold on; hold all;
% plot(mean(csensb_ccw_d,1),'linewidth',1, 'color', [0.5,0.5,1,0.5]);
% ylabel('Average signed DV (acc. bound)')
% xlabel('Time')
% title('Average over trials')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% -800 and -600 ms plotted against the same statistic computed between -500 and -300.

% ==> timepoints (like the one for Monkey F and session 1 in the physiology data)
t = linspace(-795,-45,200);

% get indices of timepoints closest to -500 and -300 ms
[~,il] = min(abs(t + 500)); %-500
[~,iu] = min(abs(t + 300)); %-300
% get indices of timepoints closest to -800 and -600 ms
[~,il2] = min(abs(t + 800)); %-800
[~,iu2] = min(abs(t + 600)); %-600

% ==> drifting prior expectation
% => time window 1
csensb_ccw_d_mu =  mean(mean(csensb_ccw_d(:,il:iu),1),2);
csensb_cw_d_mu  =  mean(mean(csensb_cw_d(:, il:iu),1),2);
% => time window 2
csensb_ccw_d_mu2 = mean(mean(csensb_ccw_d(:,il2:iu2),1),2);
csensb_cw_d_mu2  = mean(mean(csensb_cw_d(:, il2:iu2),1),2);

% ==> fixed prior expectation
% => time window 1
csensb_ccw_mu =  mean(mean(csensb_ccw(:,il:iu),1),2);
csensb_cw_mu  =  mean(mean(csensb_cw(:, il:iu),1),2);
% => time window 2
csensb_ccw_mu2 = mean(mean(csensb_ccw(:,il2:iu2),1),2);
csensb_cw_mu2  = mean(mean(csensb_cw(:, il2:iu2),1),2);

subplot(1,2,1);
hold on; hold all;
% scatter(csensb_cw_d_mu2 - csensb_ccw_d_mu2, csensb_cw_d_mu - csensb_ccw_d_mu,  90, 'ko'); %,'markerfacecolor', 'b', 'markeredgecolor', 'b')
scatter(csensb_cw_d_mu2,  csensb_cw_d_mu,     90, 'o','markerfacecolor', 'r', 'markeredgecolor', 'r')
scatter(csensb_ccw_d_mu2,  csensb_ccw_d_mu,   90, 'o','markerfacecolor', 'b', 'markeredgecolor', 'b')
axis equal;
subplot(1,2,2);
hold on; hold all;
% scatter(csensb_cw_mu2 - csensb_ccw_mu2,     csensb_cw_mu - csensb_ccw_mu, 90,  'ko'); %,'markerfacecolor', 'b', 'markeredgecolor', 'b')
scatter(csensb_cw_mu2,    csensb_cw_mu,    90,  'o','markerfacecolor', 'r', 'markeredgecolor', 'r')
scatter(csensb_ccw_mu2,    csensb_ccw_mu,  90,  'o','markerfacecolor', 'b', 'markeredgecolor', 'b')
axis equal;

fprintf('completed plotting for prior offset factor %d...\n',fc)

end

subplot(1,2,1);
hold on; hold all;
plot(-3:0.1:3,-3:0.1:3,'k--')
plot(-3:0.1:3,zeros([1,length(-3:0.1:3)]),'k--')
plot(zeros([1,length(-3:0.1:3)]),-3:0.1:3,'k--')
axis equal; axis square;
xlabel('Average (-800ms to -600ms) window')
ylabel('Average (-500ms to -300ms) window')
legend('cw \Delta DV between time windows', 'ccw \Delta DV between time windows', 'Location', 'NorthWest')
title('Diffusion sticky bound with mean drift')
drawnow;

subplot(1,2,2);
hold on; hold all;
plot(-3:0.1:3,-3:0.1:3,'k--')
plot(-3:0.1:3,zeros([1,length(-3:0.1:3)]),'k--')
plot(zeros([1,length(-3:0.1:3)]),-3:0.1:3,'k--')
axis equal; axis square;
xlabel('Average (-800ms to -600ms) window')
ylabel('Average (-500ms to -300ms) window')
legend('cw \Delta DV between time windows', 'ccw \Delta DV between time windows', 'Location', 'NorthWest')
title('Diffusion sticky bound (no drift)')
drawnow;

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

    % => if Monkey F
    if iSfln <= 13
        % ==> store averages in two time windows
        % => lo contrast       
        Fmu_cw_lo(iSfln,:) = [mean(mean(vcwlo(:,il2:iu2),1),2), mean(mean(vcwlo(:,il:iu),1),2)];
        % => hi contrast
        Fmu_cw_hi(iSfln,:) = [mean(mean(vcwhi(:,il2:iu2),1),2), mean(mean(vcwhi(:,il:iu),1),2)];        
        % => lo contrast
        Fmu_ccw_lo(iSfln,:) = [mean(mean(vccwlo(:,il2:iu2),1),2), mean(mean(vccwlo(:,il:iu),1),2)];
        % => hi contrast
        Fmu_ccw_hi(iSfln,:) = [mean(mean(vccwhi(:,il2:iu2),1),2), mean(mean(vccwhi(:,il:iu),1),2)];    
        
        % ==> store average DV trajectories        
        F_cw_lo(iSfln,:) = mean(vcwlo,1);    
        F_cw_hi(iSfln,:) = mean(vcwhi,1);  
        F_ccw_lo(iSfln,:) = mean(vccwlo,1);    
        F_ccw_hi(iSfln,:) = mean(vccwhi,1);          
        
    % => if Monkey J    
    elseif iSfln >= 14
        % ==> store averages in two time windows        
        % => lo contrast
        Jmu_cw_lo(iSfln -13,:) = [mean(mean(vcwlo(:,il2:iu2),1),2), mean(mean(vcwlo(:,il:iu),1),2)];
        % => hi contrast
        Jmu_cw_hi(iSfln -13,:) = [mean(mean(vcwhi(:,il2:iu2),1),2), mean(mean(vcwhi(:,il:iu),1),2)];        
        % => lo contrast
        Jmu_ccw_lo(iSfln -13,:) = [mean(mean(vccwlo(:,il2:iu2),1),2), mean(mean(vccwlo(:,il:iu),1),2)];
        % => hi contrast
        Jmu_ccw_hi(iSfln -13,:) = [mean(mean(vccwhi(:,il2:iu2),1),2), mean(mean(vccwhi(:,il:iu),1),2)];      
        
        % ==> store average DV trajectories        
        J_cw_lo(iSfln  -13,:) = mean(vcwlo,1);    
        J_cw_hi(iSfln  -13,:) = mean(vcwhi,1);  
        J_ccw_lo(iSfln -13,:) = mean(vccwlo,1);    
        J_ccw_hi(iSfln -13,:) = mean(vccwhi,1);    
        
    end    
end

%%
% close all; clc;

figure(); set(gcf,'Color','w'); set(gcf,'Position',[86 490 1834 472]);
subplot(1,3,1);
hold on; hold all;
% => monkey F
s1 = scatter(Fmu_cw_lo(:,1), Fmu_cw_lo(:,2),  90, 'o','markerfacecolor', [0.15,0.75,0.5], 'markeredgecolor', 'r', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s2 = scatter(Fmu_ccw_lo(:,1),Fmu_ccw_lo(:,2), 90, 'o','markerfacecolor', [0.15,0.75,0.5], 'markeredgecolor', 'b', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s3 = scatter(Fmu_cw_hi(:,1), Fmu_cw_hi(:,2),  90, 'o','markerfacecolor', [0.15,0.75,0.5], 'markeredgecolor', 'r');
s4 = scatter(Fmu_ccw_hi(:,1),Fmu_ccw_hi(:,2), 90, 'o','markerfacecolor', [0.15,0.75,0.5], 'markeredgecolor', 'b');
% => monkey J
s5 = scatter(Jmu_cw_lo(:,1), Jmu_cw_lo(:,2),  90, 'o','markerfacecolor', [1,0.5,0], 'markeredgecolor', 'r', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s6 = scatter(Jmu_ccw_lo(:,1),Jmu_ccw_lo(:,2), 90, 'o','markerfacecolor', [1,0.5,0], 'markeredgecolor', 'b', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
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

%% ==> just seprate out by monkey

figure(); set(gcf,'Color','w'); set(gcf,'Position',[178 358 1325 604]);
subplot(1,2,1);
hold on; hold all;
% => monkey F
s1 = scatter(Fmu_cw_lo(:,1), Fmu_cw_lo(:,2),  140, 'o','markerfacecolor', [0.15,0.75,0.5], 'markeredgecolor', 'r', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s2 = scatter(Fmu_ccw_lo(:,1),Fmu_ccw_lo(:,2), 140, 'o','markerfacecolor', [0.15,0.75,0.5], 'markeredgecolor', 'b', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
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
s1 = scatter(Jmu_cw_lo(:,1), Jmu_cw_lo(:,2),  140, 'o','markerfacecolor', [1,0.5,0], 'markeredgecolor', 'r', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s2 = scatter(Jmu_ccw_lo(:,1),Jmu_ccw_lo(:,2), 140, 'o','markerfacecolor', [1,0.5,0], 'markeredgecolor', 'b', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
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
xlabel('Time'); ylabel('Average Signed DV (vertical stimulus)');
ylim([-1,1]);
% => time windows
plot(repmat(-800,[20,1]),linspace(-0.8,0.8,20),'k--')
plot(repmat(-600,[20,1]),linspace(-0.8,0.8,20),'k--')
plot(repmat(-500,[20,1]),linspace(-0.8,0.8,20),'k:')
plot(repmat(-300,[20,1]),linspace(-0.8,0.8,20),'k:')
xlim([-840,-30]);

hold on; hold all;   
p3 = plot(t,mean(vcwhi,1),  'r-','linewidth',4);
p4 = plot(t,mean(vccwhi,1), 'b-','linewidth',4);    

legend([p1,p2,p3,p4],'CW prior low contrast','CCW prior low contrast', 'CW prior high contrast','CCW prior high contrast');
drawnow;
    
