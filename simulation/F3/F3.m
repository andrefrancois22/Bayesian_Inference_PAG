clear all; close all; %clc;

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
N = 1000;
% => number of simulated timepoints
tm = 200; % want total tim points to equal 200 (like actual DV fits)  % 200;
% ==> bound
bnd = 12; %; %10

% ==> dynamic range function
dynf = @(x) max([max(x, [], 2) - x(:,1), -(min(x, [], 2) - x(:,1))], [], 2);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lim = 20;
figure(10); set(gcf,'color','white'); set(gcf,'Position',[495 530 770 399]);
subplot(1,2,2);
axis equal;
hold on; hold all;
plot(-lim:0.1:lim,-lim:0.1:lim,'k--')
plot(-lim:0.1:lim,zeros([1,length(-lim:0.1:lim)]),'k--')
plot(zeros([1,length(-lim:0.1:lim)]),-lim:0.1:lim,'k--')
axis equal; axis square;
xlabel('Average (-800ms to -600ms) window')
ylabel('Average (-500ms to -300ms) window')
% legend('cw', 'ccw', 'Location', 'NorthWest')
title(['Diffusion sticky bound with mean drift '])
drawnow; 


% ==> store proportions for each dyn range split
prop_cw  = nan(2,7);
prop_ccw = nan(2,7);

% ==> standard deviation for randn
sd = 1;
% ==> mean for randn
mus = [-0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15]; %****

% ==> prior drift linear factors
fcs = 0:0.25:10; %0:0.01:1; %
fci = 1;

% ==> stimulus orientation
for fc = fcs 
    
    % ==> factor for prior offset in cw and ccw directions
    for or = 1:7

         % ==> Static prior offsets (model 1)
        pr_cw  =  25 * fc;
        pr_ccw = -25 * fc;

        % ==> randn
        rdn = randn(N,tm);
        % ==> Gaussian random walk
        sens_cw  = mus(or) + sd*rdn;
        sens_ccw = mus(or) + sd*rdn;
        % ==> cummulative sum over time
        csens_cw  = cumsum(sens_cw,2);
        csens_ccw = cumsum(sens_ccw,2);  

        % ==> linear drift case
        % ==> case with a bias that grows linearly with t (mean drift)
        pr_cw_d  =  fc*(0:(tm-1)) + pr_cw; %****
        pr_ccw_d = -fc*(0:(tm-1)) + pr_ccw; %****

        % ==> model 4
        % Linear drift starts during -800ms - -600ms time window
        pr_cw_d_m  = pr_cw_d  / tm; %(0.5*tm);
        pr_ccw_d_m = pr_ccw_d / tm; %(0.5*tm);
        % ==> csens with bound - two cases: with ccw prior or with cw prior
        csensb_cw_d =  csens_cw  + repmat(pr_cw_d_m,  [N,1]); %****
        csensb_ccw_d = csens_ccw + repmat(pr_ccw_d_m, [N,1]); %****   

%         csensb_cw_d =  csens_cw  + pr_cw;  
%         csensb_ccw_d = csens_ccw + pr_ccw;
        
        figure(10); 
        subplot(1,2,1); hold on; hold all; 
        plot(t, mean(csensb_cw_d,1),  'r-');
        plot(t, mean(csensb_ccw_d,1), 'b-');
        drawnow;
        subplot(1,2,2); hold on; hold all; 
        scatter(mean(mean(csensb_cw_d(:,il2:iu2))),mean(mean(csensb_cw_d(:,il:iu))),   'r.')
        scatter(mean(mean(csensb_ccw_d(:,il2:iu2))),mean(mean(csensb_ccw_d(:,il:iu))), 'b.')
        axis equal; drawnow;
        
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

        % ==> updates
        fprintf('completed plotting for prior offset factor %d...\n',fc)

        %==> simulating responses and dynamic range analysis
        % ==> responses for bound crossings
        % ==> cw
        cw_dv_cw_r =    (csensb_cw_d(:,end)  ==   bnd);
        cw_dv_ccw_r =  -(csensb_cw_d(:,end)  ==  -bnd);
        % ==> ccw
        ccw_dv_cw_r =   (csensb_ccw_d(:,end) ==  bnd);
        ccw_dv_ccw_r = -(csensb_ccw_d(:,end) == -bnd);
        % ==> some sanity checks
        assert(sum([cw_dv_cw_r  & cw_dv_ccw_r]) == 0)
        assert(sum([ccw_dv_cw_r & ccw_dv_ccw_r]) == 0)

        % ==> create 2AFC response vector (for cw context DVs)
        cw_r  = [cw_dv_cw_r  + cw_dv_ccw_r]; 
        % ==> for ccw context DVs
        ccw_r = [ccw_dv_cw_r + ccw_dv_ccw_r];
        
%         % ==> track DVs that never crossed a bound
%         cw_r_id  = cw_r~=0; 
%         ccw_r_id = ccw_r~=0; 
%         % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
%         % ==> handle DVs that never reached a bound...
%         % ==> coin toss (ct) for zeros (trials that never reached a bound)
%         % ==> coin toss (cw DVs)
%         cw_ct = binornd(1,0.5,length(cw_r(cw_r==0)),1);
%         cw_ct(cw_ct==0) = -1;
%         % ==> coin toss (ccw DVs)
%         ccw_ct = binornd(1,0.5,length(ccw_r(ccw_r==0)),1);
%         ccw_ct(ccw_ct==0) = -1;
%         % ==> populate trials that never reached a bound with 2AFC coin toss (guess)
%         cw_r(cw_r == 0)   = cw_ct;
%         ccw_r(ccw_r == 0) = ccw_ct;

        % ==> incongruent trials (cw context yields a ccw response, e.g. '-1' )
        dvs_i_cw = csensb_cw_d(cw_r==-1,:);
        dynr_i_cw = dynf(dvs_i_cw);
        % ==> congruent trials (cw context)
        dvs_c_cw = csensb_cw_d(cw_r== 1,:);
        dynr_c_cw = dynf(dvs_c_cw); 

        % ==> incongruent trials (ccw context yields a cw response, e.g. '-1' )
        dvs_i_ccw = csensb_ccw_d(ccw_r== 1,:);
        dynr_i_ccw = dynf(dvs_i_ccw); 
        % ==> congruent trials (ccw context yields a ccw response e.g. '-1')
        dvs_c_ccw = csensb_ccw_d(ccw_r==-1,:);
        dynr_c_ccw = dynf(dvs_c_ccw); 

        % ==> cw context, and stim orientation dynr median
        med_cw  = median(dynf(csensb_cw_d));
        % ==> ccw context, and stim orientation dynr median
        med_ccw = median(dynf(csensb_ccw_d));

        % ==> proportions => high dynamic range DVs
        prop_cw(1,or)  = sum(cw_r(dynf(csensb_cw_d) > med_cw) == 1)   / (sum(cw_r(dynf(csensb_cw_d) > med_cw) == 1)  + sum(cw_r(dynf(csensb_cw_d) > med_cw) == -1));
        prop_ccw(1,or) = sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == 1)  / (sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == 1) + sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == -1));

        % ==> proportions => low dynamic range DVs
        prop_cw(2,or)  = sum(cw_r(dynf(csensb_cw_d) <= med_cw) == 1)   / (sum(cw_r(dynf(csensb_cw_d) <= med_cw) == 1)  + sum(cw_r(dynf(csensb_cw_d) <= med_cw) == -1));
        prop_ccw(2,or) = sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == 1)  / (sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == 1) + sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == -1));

    end
    
    % ==> plot simulated choice proportions by dynamic range split
    figure(3); set(gcf,'color','white'); set(gcf, 'Position',[273 218 1150 719]);
    subplot(6,7,fci)
    hold on; hold all; 
    plot([mus],[prop_cw(1,:)],'r*-'); 
    plot([mus],[prop_ccw(1,:)],'b*-');  
    plot([mus],[prop_cw(2,:)],'r.:'); 
    plot([mus],[prop_ccw(2,:)],'b.:'); 
    xlim([-max(mus),max(mus)]); ylim([0,1]);
    drawnow;
    
    fci = fci + 1;
    fprintf('finished stimulus orientation %d...\n',mus(or))
end


figure(); set(gcf,'color','white'); set(gcf,'Position',[273 211 971 726])
subplot(1,2,1);
hold on; hold all;
% ==> cw context
plot(t,mean(dvs_i_cw,1)','r--'); 
plot(t,mean(dvs_c_cw,1)','r-')
scatter(t(1),mean(dvs_i_cw(:,1),1)', 40, 'ro'); 
scatter(t(1),mean(dvs_c_cw(:,1),1)', 40, 'ro','filled');
% ==> ccw context
plot(t,mean(dvs_i_ccw,1)','b--'); 
plot(t,mean(dvs_c_ccw,1)','b-'); 
scatter(t(1),mean(dvs_i_ccw(:,1),1)', 40, 'bo'); 
scatter(t(1),mean(dvs_c_ccw(:,1),1)', 40, 'bo','filled');
plot(t,zeros(200,1),'k--')
legend('cw incongruent','cw congruent', ...
       'cw incongruent DV start','cw congruent DV start', ...
       'ccw incongruent','ccw congruent', ...
       'ccw incongruent DV start','ccw congruent DV start', ...
       'Fontsize',6)
xlim([t(1)-50,t(end)])   
ylim([-20,20]);
xlabel('Time');
ylabel('Simulated DV (average)');
title(['Simulated DV (average) stimulus: ',num2str(or)])


ax = subplot(1,2,2);
hold  on; hold all;
errorbar([1,2],[mean(dynr_i_cw), mean(dynr_c_cw)], [std(dynr_i_cw), std(dynr_c_cw)],'r.-')
errorbar([3,4],[mean(dynr_i_ccw),mean(dynr_c_ccw)],[std(dynr_i_ccw),std(dynr_c_ccw)],'b.-')
xlim([0,5]);
xlabel('Context and choice congruence');
ylabel('Dynamic range (average)');
ax.XTick = [1,2,3,4];
ax.XTickLabel = {'cw (incongruent)', 'cw (congruent)', 'ccw (incongruent)', 'ccw (congruent)'};
ax.XTickLabelRotation = 45;
