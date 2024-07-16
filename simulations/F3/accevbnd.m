function [db,dp, prop_cw, prop_ccw, dvs_c_cw, dvs_i_cw, dvs_c_ccw, dvs_i_ccw] = accevbnd(N, tm, bnd, sd, mus, fcs, ofs, M_FLAG)

% ==> delta bias
db = nan(1,length(fcs));
% ==> delta perceptual uncertainty
dp = nan(1,length(fcs));

% ==> store proportions for each dyn range split
prop_cw  = nan(length(fcs),3,7);
prop_ccw = nan(length(fcs),3,7);

% ==> dynamic range function
dynf = @(x) max([max(x, [], 2) - x(:,1), -(min(x, [], 2) - x(:,1))], [], 2);

fci = 1;

for fc = fcs
        
    % ==> stimulus orientation
    for or = 1:7

        % ==> factor for prior offset in cw and ccw directions
        % ==> prior offsets
        pr_cw  =  ofs * fc; 
        pr_ccw = -ofs * fc; 
        
        pr_cw_v  = repmat(pr_cw,  [tm,1]);  pr_cw_v(1)  =  ofs;
        pr_ccw_v = repmat(pr_ccw, [tm,1]);  pr_ccw_v(1) = -ofs;
        
        % ==> randn
        rdn = randn(N,tm);
        % ==> Gaussian random walk
        sens_cw  = mus(or) + sd*rdn; %(momentary sens evidence)
        sens_ccw = mus(or) + sd*rdn;        
        
        % ==> FIX HERE..
        % ==> cummulative sum over time (momentary evidence - )
        csensb_cw_d  = cumsum(sens_cw  + repmat(pr_cw_v', [N,1]), 2);
        csensb_ccw_d = cumsum(sens_ccw + repmat(pr_ccw_v', [N,1]),2);          

        % ==> once a bnd has been reached at time t for a given trial, set the remaining DV values
        % starting from t+1 to the bound for that trial
        % ==> cw cases
        csensb_cw_d = bndf(csensb_cw_d,bnd,N, M_FLAG);
        % ==> ccw context cases
        csensb_ccw_d = bndf(csensb_ccw_d,bnd,N, M_FLAG);

        if strcmp(M_FLAG,'BOUND')
            %==> simulating responses and dynamic range analysis
            % ==> responses for bound crossings
            % ==> cw
            cw_dv_cw_r =    (csensb_cw_d(:,end)  ==   bnd);
            cw_dv_ccw_r =  -(csensb_cw_d(:,end)  ==  -bnd);
            % ==> ccw
            ccw_dv_cw_r =   (csensb_ccw_d(:,end) ==  bnd);
            ccw_dv_ccw_r = -(csensb_ccw_d(:,end) == -bnd);
            
        elseif strcmp(M_FLAG,'NO_BOUND')   
            % ==> no bound - responses are sign at the end of DV
            % ==> cw
            cw_dv_cw_r =    (sign(csensb_cw_d(:,end))  ==   1);
            cw_dv_ccw_r =  -(sign(csensb_cw_d(:,end))  ==  -1);
            % ==> ccw
            ccw_dv_cw_r =   (sign(csensb_ccw_d(:,end)) ==  1);
            ccw_dv_ccw_r = -(sign(csensb_ccw_d(:,end)) == -1);                        
        end              
        
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
        % ==> congruent trials (cw context)
        dvs_c_cw = csensb_cw_d(cw_r== 1,:);
 
        % ==> incongruent trials (ccw context yields a cw response, e.g. '-1' )
        dvs_i_ccw = csensb_ccw_d(ccw_r== 1,:);
        % ==> congruent trials (ccw context yields a ccw response e.g. '-1')
        dvs_c_ccw = csensb_ccw_d(ccw_r==-1,:);

        % ==> cw context, and stim orientation dynr median
        med_cw  = median(dynf(csensb_cw_d));
        % ==> ccw context, and stim orientation dynr median
        med_ccw = median(dynf(csensb_ccw_d));

        % ==> proportions => high dynamic range DVs
        prop_cw(fci,2,or)  = sum(cw_r(dynf(csensb_cw_d) > med_cw) == 1)     / (sum(cw_r(dynf(csensb_cw_d) > med_cw) == 1)    + sum(cw_r(dynf(csensb_cw_d) > med_cw) == -1));
        prop_ccw(fci,2,or) = sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == 1)  / (sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == 1) + sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == -1));

        % ==> proportions => low dynamic range DVs
        prop_cw(fci,1,or)  = sum(cw_r(dynf(csensb_cw_d) <= med_cw) == 1)     / (sum(cw_r(dynf(csensb_cw_d) <= med_cw) == 1)    + sum(cw_r(dynf(csensb_cw_d) <= med_cw) == -1));
        prop_ccw(fci,1,or) = sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == 1)  / (sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == 1) + sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == -1));

        % ==> proportions without split
        prop_cw(fci,3,or)  = sum(cw_r  == 1) / (sum(cw_r  == 1) + sum(cw_r  == -1));
        prop_ccw(fci,3,or) = sum(ccw_r == 1) / (sum(ccw_r == 1) + sum(ccw_r == -1));

    end
    
%     figure(2);
%     plot(1:tm,mean(csensb_cw_d(cw_r == 1,:)),'r-')
%     hold on; hold all;
%     plot(1:tm,mean(csensb_ccw_d(ccw_r == -1,:)),'b-')    
%     keyboard
        
    % ==> delta bias
    db(fci) = (prop_cw(fci,2,4) - prop_ccw(fci,2,4)) - (prop_cw(fci,1,4) - prop_ccw(fci,1,4));
    % ==> delta perceptual uncertainty (approx)
    dp(fci) = (prop_cw(fci,1,5) - prop_cw(fci,1,3)) - (prop_cw(fci,2,5) - prop_cw(fci,2,3));
        
    % ==> update
    fprintf('Simulated PFs for drift rate and offset parameter %d of %d...\n', fci, length(fcs));
    % ==> increment
    fci = fci + 1;    
end
