function [db,dp] = accevbnd(N, tm, bnd, sd, mus, fcs)

% ==> delta bias
db = nan(1,length(fcs));
% ==> delta perceptual uncertainty
dp = nan(1,length(fcs));

% ==> store proportions for each dyn range split
prop_cw  = nan(3,7);
prop_ccw = nan(3,7);

% ==> dynamic range function
dynf = @(x) max([max(x, [], 2) - x(:,1), -(min(x, [], 2) - x(:,1))], [], 2);

fci = 1;
% ==> stimulus orientation
for fc = fcs 
    
    % ==> factor for prior offset in cw and ccw directions
    for or = 1:7

         % ==> prior offsets
        pr_cw  =  15 * fc; %25
        pr_ccw = -15 * fc; %-25
        % ==> linear drift case
        % ==> case with a bias that grows linearly with t (mean drift)
        pr_cw_d  =  fc*(0:(tm-1)) + pr_cw;  
        pr_ccw_d = -fc*(0:(tm-1)) + pr_ccw;

        % ==> randn
        rdn = randn(N,tm);
        % ==> Gaussian random walk
        sens_cw  = mus(or) + sd*rdn;
        sens_ccw = mus(or) + sd*rdn;
        % ==> cummulative sum over time
        csens_cw  = cumsum(sens_cw,2);
        csens_ccw = cumsum(sens_ccw,2);          
        
        % Linear drift starts during -800ms - -600ms time window
        pr_cw_d_m  = pr_cw_d  / tm;
        pr_ccw_d_m = pr_ccw_d / tm; 
        % ==> csens with bound - two cases: with ccw prior or with cw prior
        csensb_cw_d =  csens_cw  + repmat(pr_cw_d_m,  [N,1]); 
        csensb_ccw_d = csens_ccw + repmat(pr_ccw_d_m, [N,1]);    

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
        csensb_cw_d = bndf(csensb_cw_d,bnd,N);
        % ==> ccw context cases
        csensb_ccw_d = bndf(csensb_ccw_d,bnd,N);

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
        
        % ==> handle DVs that never reached a bound...
        % ==> use sign of DV
        cw_r(cw_r == 0)   = sign(csensb_cw_d(cw_r==0,end));
        ccw_r(ccw_r == 0) = sign(csensb_ccw_d(ccw_r==0,end));

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
        prop_cw(1,or)  = sum(cw_r(dynf(csensb_cw_d) > med_cw) == 1)     / (sum(cw_r(dynf(csensb_cw_d) > med_cw) == 1)    + sum(cw_r(dynf(csensb_cw_d) > med_cw) == -1));
        prop_ccw(1,or) = sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == 1)  / (sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == 1) + sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == -1));

        % ==> proportions => low dynamic range DVs
        prop_cw(2,or)  = sum(cw_r(dynf(csensb_cw_d) <= med_cw) == 1)     / (sum(cw_r(dynf(csensb_cw_d) <= med_cw) == 1)    + sum(cw_r(dynf(csensb_cw_d) <= med_cw) == -1));
        prop_ccw(2,or) = sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == 1)  / (sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == 1) + sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == -1));

        % ==> proportions without split
        prop_cw(3,or)  = sum(cw_r  == 1) / (sum(cw_r  == 1) + sum(cw_r  == -1));
        prop_ccw(3,or) = sum(ccw_r == 1) / (sum(ccw_r == 1) + sum(ccw_r == -1));

    end
    
    % ==> delta bias
    db(fci) = (prop_cw(2,4) - prop_ccw(2,4)) - (prop_cw(1,4) - prop_ccw(1,4));
    % ==> delta perceptual uncertainty (approx)
    dp(fci) = (prop_cw(1,5) - prop_cw(1,3)) - (prop_cw(2,5) - prop_cw(2,3));
    
    % ==> plot simulated choice proportions by dynamic range split
    figure(3); set(gcf,'color','white'); set(gcf, 'Position',[273 218 1150 719]);
    subplot(6,7,fci)
    hold on; hold all; 
    plot([mus],[prop_cw(3,:)],'b.-', 'linewidth', 1.5); 
    plot([mus],[prop_ccw(3,:)],'r.-', 'linewidth', 1.5); 
    xlabel('Orientation');
    ylabel('p(cw)')    
    %xlim([-max(mus),max(mus)]); 
    ylim([0,1]);
    axis square;
    drawnow;      
    % ==> plot simulated choice proportions by dynamic range split
    figure(4); set(gcf,'color','white'); set(gcf, 'Position',[273 218 1150 719]);
    subplot(6,7,fci)
    hold on; hold all; 
    plot([mus],[prop_cw(1,:)],'bx:', 'linewidth', 1); 
    plot([mus],[prop_ccw(1,:)],'rx:','linewidth', 1);  
    plot([mus],[prop_cw(2,:)],'b.-', 'linewidth', 1.5); 
    plot([mus],[prop_ccw(2,:)],'r.-', 'linewidth', 1.5); 
    xlabel('Orientation');
    ylabel('p(cw)')
    %xlim([-max(mus),max(mus)]); 
    ylim([0,1]);
    axis square;
    drawnow;      
    % ==> update
    fprintf('imulated PFs for drift rate and offset parameter %d of %d...\n', fci, length(fcs));
    fci = fci + 1;    
end