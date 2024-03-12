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

mdl = 'model 3';

% ==> figure 
figure(1); set(gcf,'Color','w'); set(gcf, 'Position',[71 510 1850 440])

% ==> factor for prior offset in cw and ccw directions
for fc = 1:0.5:10 %1:0.1:5

     % ==> Static prior offsets (model 1)
    pr_cw  =  0.5 * fc;
    pr_ccw = -0.5 * fc;

    % ==> Gaussian random walk
    sens  = randn(N,tm);
    % ==> cummulative sum over time
    csens = cumsum(sens,2);
    % ==> csens with bound -> two cases: with ccw prior or with cw prior offsets
    csensb_cw =  csens + pr_cw;
    csensb_ccw = csens + pr_ccw;
    
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

    % ==> model with prior offset and linear drift (model 3)
    % ==> case with a bias that grows linearly with t (mean drift)
    pr_cw_d  =  pr_cw  +fc*(0:(tm-1));
    pr_ccw_d =  pr_ccw -fc*(0:(tm-1));

    if strcmp(mdl, 'model 2')
        % ==> model 3
        pr_cw_d_m  = [repmat(pr_cw, [1,iu2]),  pr_cw_d(1:(tm - iu2))]  /(0.5*tm);
        pr_ccw_d_m = [repmat(pr_ccw, [1,iu2]), pr_ccw_d(1:(tm - iu2))] /(0.5*tm);  
        % ==> csens with bound - two cases: with ccw prior or with cw prior
        csensb_cw_d =  csens + repmat(pr_cw_d_m,  [N,1]); %repmat(pr_cw_d, [N,1]);
        csensb_ccw_d = csens + repmat(pr_ccw_d_m, [N,1]); %repmat(pr_ccw_d,[N,1]);    
    end

    if strcmp(mdl, 'model 3')
        % ==> model 2
        pr_cw_d_m  = pr_cw_d  / (0.5*tm);
        pr_ccw_d_m = pr_ccw_d / (0.5*tm);
        % ==> csens with bound - two cases: with ccw prior or with cw prior
        csensb_cw_d =  csens + repmat(pr_cw_d_m,  [N,1]); %repmat(pr_cw_d, [N,1]);
        csensb_ccw_d = csens + repmat(pr_ccw_d_m, [N,1]); %repmat(pr_ccw_d,[N,1]);    
    end
    
    % ==> plot
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
    subplot(1,3,3);
    hold on; hold all;
    scatter(csensb_cw_d_mu2,  csensb_cw_d_mu,     90, 'o','markerfacecolor', 'r', 'markeredgecolor', 'r')
    scatter(csensb_ccw_d_mu2,  csensb_ccw_d_mu,   90, 'o','markerfacecolor', 'b', 'markeredgecolor', 'b')
    axis equal;
    
    % ==> plot simulated results without drift
    subplot(1,3,2);
    hold on; hold all;
    scatter(csensb_cw_mu2,    csensb_cw_mu,    90,  'o','markerfacecolor', 'r', 'markeredgecolor', 'r')
    scatter(csensb_ccw_mu2,    csensb_ccw_mu,  90,  'o','markerfacecolor', 'b', 'markeredgecolor', 'b')
    axis equal;
    
    % ==> updates
    fprintf('completed plotting for prior offset factor %d...\n',fc)

end

lim = fc/2;
% ==> plot axes and lines
subplot(1,3,3);
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

% ==> plot axes and lines
subplot(1,3,2);
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
