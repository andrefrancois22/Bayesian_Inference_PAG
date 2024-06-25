function [nll] = calcnllf(iS,cr,CTs,PF,propccw,propcw)
% ==> compute negative-log-likelihood of simulation dynamic range PF curve estimates

% ==> initialize vector of log-likelihoods for 4 simulated PFs curves
% <dynamic range split>_<context>_<stimulus orientation>
lls = nan(2,2,7);
for sp = 1:2
%     % ==> compute lls (7 values for each of the 4 conditions - 2 x 2 with dyn range split, and context)
%     lls(sp,1,:) = log(max(1e-300, binopdf( CTs{iS,cr,sp,2}(1,:), CTs{iS,cr,sp,1}(1,:) + CTs{iS,cr,sp,2}(1,:), propccw(sp,:))));    
%     lls(sp,2,:) = log(max(1e-300, binopdf( CTs{iS,cr,sp,2}(2,:), CTs{iS,cr,sp,1}(2,:) + CTs{iS,cr,sp,2}(2,:), propcw(sp,:))));
%     
    % ==> try SSE
    lls(sp,1,:) = (PF{iS,cr,sp}(2,:) - propccw(sp,:)).^2;    
    lls(sp,2,:) = (PF{iS,cr,sp}(1,:) - propcw(sp,:)).^2;    

end
% ==> compute sum (minimizing this quantity amounts to maximizing the ll)
% nll = -sum(lls(:));
nll = sum(lls(:));
end