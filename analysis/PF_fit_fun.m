function [pf_fit, pf_params] = PF_fit_fun(nccw_inc, nccw_cng, ncw_cng, ncw_inc, or, startVec_M1, LB_M1, UB_M1, options)
    % function for fitting PF curve for each context using 2AFC counts.

    % ==> counts
    % ==> context 1 ccw response (incongruent) & ==> context 2 ccw response (congruent)
    nccw_dynr_lo = [nccw_inc; nccw_cng]; 
    % ==> context 1 cw response (congruent)    & ==> context 2 cw response (incongruent)
    ncw_dynr_lo  = [ncw_cng;  ncw_inc]; 
    
    % ==> pass through PF function fit
    obFun    = @(paramVec) giveNLL(paramVec, or, nccw_dynr_lo, ncw_dynr_lo, 'M1');    
    % ==> run 
    pf_params = fmincon(obFun, startVec_M1, [], [], [], [], LB_M1, UB_M1, [], options);
    % ==> return pf
    [~, pf_fit] = giveNLL(pf_params, or, nccw_dynr_lo, ncw_dynr_lo, 'M1');

    
end

