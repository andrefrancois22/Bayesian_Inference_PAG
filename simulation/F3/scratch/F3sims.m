clear all; close all; clc;

% ==> parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> timepoints (like the one for Monkey F and session 1 in the physiology data)
t = linspace(-795,-45,200);
% get indices of timepoints closest to -500 and -300 ms
[~,il] = min(abs(t + 500)); %-500
[~,iu] = min(abs(t + 300)); %-300
% get indices of timepoints closest to -800 and -600 ms
[~,il2] = min(abs(t + 800)); %-800
[~,iu2] = min(abs(t + 600)); %-600

% => number of simulated trials
N = 500;
% => number of simulated timepoints
tm = 200; % want total tim points to equal 200 (like actual DV fits)  % 200;
% ==> bound
bnd = 10;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> choose model
% mdl = 'model2';
% mdl = 'model3';
mdl = 'model4';

% ==> standard deviation for randn
sd = 1;
% ==> mean for randn
mus = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]; %****
% ==> figure 
figure(1); set(gcf,'Color','w'); set(gcf, 'Position',[71 510 1850 440])

% ==> factor for prior offset in cw and ccw directions
for fc = 0:0.5:10 %1:0.1:5

     % ==> Static prior offsets (model 1)
    pr_cw  =  0.5 * fc;
    pr_ccw = -0.5 * fc;

    % ==> randn
    rdn = randn(N,tm);
    % ==> Gaussian random walk
    sens_cw  =  mus(4) + sd*rdn;
    sens_ccw = -mus(4) + sd*rdn;
    % ==> cummulative sum over time
    csens_cw  = cumsum(sens_cw,2);
    csens_ccw = cumsum(sens_ccw,2);
    
    % ==> csens with bound -> two cases: with ccw prior or with cw prior offsets
    csensb_cw =  csens_cw  + pr_cw; %****
    csensb_ccw = csens_ccw + pr_ccw; %****
    
    % ==> for model 1
    % ==> set DV values for timepoints over trials that reach or exceed a bound
    % to the bound value (stickiness)
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
        bidx = find(csensb_cw(i,:) >=  bnd, 1, 'first');
        % => set to bound
        if ~isempty(bidx)
            csensb_cw(i,bidx:end) = bnd;
        end
    end
    for i = 1:N
        bidx = find(csensb_cw(i,:) <=  -bnd);
        bidx = min(bidx);
        % => set to bound
        if ~isempty(bidx)
            csensb_cw(i,bidx:end) = -bnd;
        end
    end
    % ==> ccw context cases
    for i = 1:N
        bidx = find(csensb_ccw(i,:) >=  bnd, 1, 'first');
        % => set to bound
        if ~isempty(bidx)
            csensb_ccw(i,bidx:end) = bnd;
        end
    end
    for i = 1:N
        bidx = find(csensb_ccw(i,:) <=  -bnd);
        bidx = min(bidx);
        % => set to bound
        if ~isempty(bidx)
            csensb_ccw(i,bidx:end) = -bnd;
        end
    end

    % ==> linear drift case
    % ==> case with a bias that grows linearly with t (mean drift)
    pr_cw_d  =  + fc*(0:(tm-1));
    pr_ccw_d =  - fc*(0:(tm-1));

    if strcmp(mdl, 'model2')
        % ==> model 2
        % 0 offset during -800ms - -600ms time window, linear drift starts
        % at -600ms.
        pr_cw_d_m  = [zeros([1,iu2]),  pr_cw_d(1:(tm - iu2))]   / (0.5*tm);
        pr_ccw_d_m = [zeros([1,iu2]),  pr_ccw_d(1:(tm - iu2))]  / (0.5*tm);  
        % ==> csens with bound - two cases: with ccw prior or with cw prior
        csensb_cw_d =  csens_cw  + repmat(pr_cw_d_m,  [N,1]);
        csensb_ccw_d = csens_ccw + repmat(pr_ccw_d_m, [N,1]);    
    end
    if strcmp(mdl, 'model3')
        % ==> model 3
        % no drift (constant offset) during -800ms - -600ms time window, linear
        % drift starts at -600ms.
        pr_cw_d_m  = [repmat(pr_cw_d(iu2),  [1,iu2-1]),  pr_cw_d(iu2:end)]   / (0.5*tm);
        pr_ccw_d_m = [repmat(pr_ccw_d(iu2), [1,iu2-1]),  pr_ccw_d(iu2:end)]  / (0.5*tm);  
        % ==> csens with bound - two cases: with ccw prior or with cw prior
        csensb_cw_d =  csens_cw  + repmat(pr_cw_d_m,  [N,1]); 
        csensb_ccw_d = csens_ccw + repmat(pr_ccw_d_m, [N,1]);   
    end
    if strcmp(mdl, 'model4')
        % ==> model 4
        % Linear drift starts during -800ms - -600ms time window
        pr_cw_d_m  = pr_cw_d  / (0.5*tm);
        pr_ccw_d_m = pr_ccw_d / (0.5*tm);
        % ==> csens with bound - two cases: with ccw prior or with cw prior
        csensb_cw_d =  csens_cw  + repmat(pr_cw_d_m,  [N,1]); %****
        csensb_ccw_d = csens_ccw + repmat(pr_ccw_d_m, [N,1]); %****   
    end    
    
    % ==> plot the mean drift pattern
    % before adding cummulative gaussian noise trial trajectories)
    subplot(1,3,1);
    plot(t,pr_cw_d_m)
    hold on; hold all
    plot(t,pr_ccw_d_m)
    if fc == 10
        plot(repmat(t(iu), [1,length(-2*fc:0.1:2*fc)]),  -2*fc:0.1:2*fc,'r.-')
        plot(repmat(t(il), [1,length(-2*fc:0.1:2*fc)]),  -2*fc:0.1:2*fc,'r.-')
        plot(repmat(t(iu2),[1,length(-2*fc:0.1:2*fc)]),  -2*fc:0.1:2*fc,'k.-')
        plot(repmat(t(il2),[1,length(-2*fc:0.1:2*fc)]),  -2*fc:0.1:2*fc,'k.-')    
    end
    drawnow;
        

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
        bidx = find(csensb_cw_d(i,:) >=  bnd, 1, 'first');
        % => set to bound
        if ~isempty(bidx)
            csensb_cw_d(i,bidx:end) = bnd;
        end
    end
    for i = 1:N
        bidx = find(csensb_cw_d(i,:) <=  -bnd);
        bidx = min(bidx);
        % => set to bound
        if ~isempty(bidx)
            csensb_cw_d(i,bidx:end) = -bnd;
        end
    end
    % ==> ccw context cases
    for i = 1:N
        bidx = find(csensb_ccw_d(i,:) >=  bnd, 1, 'first');
        % => set to bound
        if ~isempty(bidx)
            csensb_ccw_d(i,bidx:end) = bnd;
        end
    end
    for i = 1:N
        bidx = find(csensb_ccw_d(i,:) <=  -bnd);
        bidx = min(bidx);
        % => set to bound
        if ~isempty(bidx)
            csensb_ccw_d(i,bidx:end) = -bnd;
        end
    end

    % ==> drifting prior expectation
    % => average over trials, and within time window 1
    csensb_ccw_d_mu =  mean(mean(csensb_ccw_d(:,il:iu),1),2);
    csensb_cw_d_mu  =  mean(mean(csensb_cw_d(:, il:iu),1),2);
    % => average over trials, and within time window 2
    csensb_ccw_d_mu2 = mean(mean(csensb_ccw_d(:,il2:iu2),1),2);
    csensb_cw_d_mu2  = mean(mean(csensb_cw_d(:, il2:iu2),1),2);

    % ==> fixed prior expectation
    % => => average over trials, and within time window 1
    csensb_ccw_mu =  mean(mean(csensb_ccw(:,il:iu),1),2);
    csensb_cw_mu  =  mean(mean(csensb_cw(:, il:iu),1),2);
    % => => average over trials, and within time window 2
    csensb_ccw_mu2 = mean(mean(csensb_ccw(:,il2:iu2),1),2);
    csensb_cw_mu2  = mean(mean(csensb_cw(:, il2:iu2),1),2);

    % ==> plot simulated results with drift
    lim = fc/2;
    subplot(1,3,3);
    hold on; hold all;
    scatter(csensb_cw_d_mu2,  csensb_cw_d_mu,     90, 'o','markerfacecolor', 'r', 'markeredgecolor', 'r')
    scatter(csensb_ccw_d_mu2,  csensb_ccw_d_mu,   90, 'o','markerfacecolor', 'b', 'markeredgecolor', 'b')
  
    
    % ==> plot simulated results without drift
    subplot(1,3,2);
    hold on; hold all;
    scatter(csensb_cw_mu2,    csensb_cw_mu,    90,  'o','markerfacecolor', 'r', 'markeredgecolor', 'r')
    scatter(csensb_ccw_mu2,    csensb_ccw_mu,  90,  'o','markerfacecolor', 'b', 'markeredgecolor', 'b')

    
    % ==> updates
    fprintf('completed plotting for prior offset factor %d...\n',fc)

end

subplot(1,3,3);
axis equal;
hold on; hold all;
plot(-lim:0.1:lim,-lim:0.1:lim,'k--')
plot(-lim:0.1:lim,zeros([1,length(-lim:0.1:lim)]),'k--')
plot(zeros([1,length(-lim:0.1:lim)]),-lim:0.1:lim,'k--')
axis equal; axis square;
xlabel('Average (-800ms to -600ms) window')
ylabel('Average (-500ms to -300ms) window')
legend('cw', 'ccw', 'Location', 'NorthWest')
title(['Diffusion sticky bound with mean drift ', mdl])
drawnow;  

subplot(1,3,2);
axis equal;
hold on; hold all;
plot(-lim:0.1:lim,-lim:0.1:lim,'k--')
plot(-lim:0.1:lim,zeros([1,length(-lim:0.1:lim)]),'k--')
plot(zeros([1,length(-lim:0.1:lim)]),-lim:0.1:lim,'k--')
axis equal; axis square;
xlabel('Average (-800ms to -600ms) window')
ylabel('Average (-500ms to -300ms) window')
legend('cw', 'ccw', 'Location', 'NorthWest')
title('Diffusion sticky bound (no drift - model 1)')
drawnow;    


%% simulating responses and dynamic range analysis

% ==> responses
cw_dv_cw_r =    (csensb_cw(:,end)  ==   bnd);
cw_dv_ccw_r =  -(csensb_cw(:,end)  ==  -bnd);

ccw_dv_cw_r =   (csensb_ccw(:,end) ==  bnd);
ccw_dv_ccw_r = -(csensb_ccw(:,end) == -bnd);

% ==> some sanity checks
assert(sum([cw_dv_cw_r  & cw_dv_ccw_r]) == 0)
assert(sum([ccw_dv_cw_r & ccw_dv_ccw_r]) == 0)

% ==> create 2AFC response vector (for cw context DVs)
cw_r  = [cw_dv_cw_r  + cw_dv_ccw_r];
% ==> for ccw context DVs
ccw_r = [ccw_dv_cw_r + ccw_dv_ccw_r];

% ==> coin toss (ct) for zeros (trials that never reached a bound)
% ==> coin toss (cw DVs)
cw_ct = binornd(1,0.5,length(cw_r(cw_r==0)),1);
cw_ct(cw_ct==0) = -1;
% ==> coin toss (ccw DVs)
ccw_ct = binornd(1,0.5,length(ccw_r(ccw_r==0)),1);
ccw_ct(ccw_ct==0) = -1;

% ==> populate trials that never reached a bound with 2AFC coin toss (guess)
cw_r(cw_r == 0)   = cw_ct;
ccw_r(ccw_r == 0) = ccw_ct;

% ==> incongruent trials (cw context yields a ccw response, e.g. '-1' )
dvs_i_cw = csensb_cw(cw_r==-1,:);
dynr_i_cw = max([max(dvs_i_cw, [], 2) - dvs_i_cw(:,1), -(min(dvs_i_cw, [], 2) - dvs_i_cw(:,1))], [], 2);
% ==> congruent trials (cw context)
dvs_c_cw = csensb_cw(cw_r== 1,:);
dynr_c_cw = max([max(dvs_c_cw, [], 2) - dvs_c_cw(:,1), -(min(dvs_c_cw, [], 2) - dvs_c_cw(:,1))], [], 2);

% ==> incongruent trials (ccw context yields a cw response, e.g. '-1' )
dvs_i_ccw = csensb_ccw(ccw_r== 1,:);
dynr_i_ccw = max([max(dvs_i_ccw, [], 2) - dvs_i_ccw(:,1), -(min(dvs_i_ccw, [], 2) - dvs_i_ccw(:,1))], [], 2);
% ==> congruent trials (ccw context yields a ccw response e.g. '-1')
dvs_c_ccw = csensb_ccw(ccw_r==-1,:);
dynr_c_ccw = max([max(dvs_c_ccw, [], 2) - dvs_c_ccw(:,1), -(min(dvs_c_ccw, [], 2) - dvs_c_ccw(:,1))], [], 2);


figure(); set(gcf,'color','white'); set(gcf, 'Position', [382 304 1171 646]);
subplot(1,2,1);
h1 = histogram(dynr_i_cw,'facecolor','r','BinWidth',1); hold on; hold all;
h2 = histogram(dynr_c_cw,'facecolor','b','BinWidth',1); 
xlabel('dynamic range')
ylabel('Frequency')
title('Dynamic range and choice congruency (cw context)')
legend([h1,h2],'incongruent trials','congruent trials')
subplot(1,2,2);
h3 = histogram(dynr_i_ccw,'facecolor','r','BinWidth',1); hold on; hold all;
h4 = histogram(dynr_c_ccw,'facecolor','b','BinWidth',1); 
xlabel('dynamic range')
ylabel('Frequency')
title('Dynamic range and choice congruency (ccw context)')
legend([h3,h4],'incongruent trials','congruent trials')
