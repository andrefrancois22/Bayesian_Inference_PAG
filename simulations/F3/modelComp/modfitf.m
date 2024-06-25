function [nll, propcw, propccw] = modfitf(iS, cr, CTs, PF, N, tm, sd, params, mod_type)
    % ==> fits accumulation to bound simple drift diffusion model
    % ==> returns simulated choice proportions for simulated dynamic range split data 
    % and calculates the negative-log-likelihood of the results 
    % compared to actual dynamic range split choice proportions
    
    % ==> stimulus orientation offset interval (assumption is that they are evenly spaced)
    intv = params(1);
    % ==> overall bias (horizontal offset of all four curves)
    bs = params(2);    
    % ==> drift rate
    fc = params(3);
    % ==> bound
    bnd = params(5);
    
    if strcmp(mod_type,'m0')     % ==> no initial offset 
        % ==> initial offset
        ofs = 0;
    elseif strcmp(mod_type,'m1') % ==> initial offset is correlated to drift rate
        % ==> initial ofset
        ofs = params(4)*fc;
    elseif strcmp(mod_type,'m2') % ==> arbitrary initial offset
        % ==> initial ofset
        ofs = params(4);        
    end
    
    % ==> orientation intervals
    mus = [-3*intv, -2*intv, -intv, 0, intv, 2*intv, 3*intv] + bs;    
    
    % ==> run accumulation of evidence to bound drift diffusion (forward) model
    % => compute simulated dynamic range split, and PFs
    [propcw, propccw] = accevbndf(N, tm, bnd, sd, mus, fc, ofs);
    
    % ==> compute nll
    [nll] = calcnllf(iS,cr,CTs, PF, propccw, propcw);
end

