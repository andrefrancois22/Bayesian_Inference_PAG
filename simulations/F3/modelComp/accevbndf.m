function [propcw, propccw] = accevbndf(N, tm, bnd, sd, mus, fc, ofs)

M_FLAG = 'BOUND';
% ==> store proportions for each dyn range split
propcw  = nan(2,7);
propccw = nan(2,7);

% ==> dynamic range function
dynf = @(x) max([max(x, [], 2) - x(:,1), -(min(x, [], 2) - x(:,1))], [], 2);
    
% ==> factor for prior offset in cw and ccw directions
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

    % ==> cw context, and stim orientation dynr median
    med_cw  = median(dynf(csensb_cw_d));
    % ==> ccw context, and stim orientation dynr median
    med_ccw = median(dynf(csensb_ccw_d));

    % ==> proportions => high dynamic range DVs
    propcw(2,or)  = sum(cw_r(dynf(csensb_cw_d) > med_cw) == 1)     / (sum(cw_r(dynf(csensb_cw_d) > med_cw) == 1)    + sum(cw_r(dynf(csensb_cw_d) > med_cw) == -1));
    propccw(2,or) = sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == 1)  / (sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == 1) + sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == -1));

    % ==> proportions => low dynamic range DVs
    propcw(1,or)  = sum(cw_r(dynf(csensb_cw_d) <= med_cw) == 1)     / (sum(cw_r(dynf(csensb_cw_d) <= med_cw) == 1)    + sum(cw_r(dynf(csensb_cw_d) <= med_cw) == -1));
    propccw(1,or) = sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == 1)  / (sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == 1) + sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == -1));
end
