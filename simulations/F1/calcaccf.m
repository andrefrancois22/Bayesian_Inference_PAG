function [acc_map, acc_mle] = calcaccf(s0,vidx,map,mle,c)
% This calculates the accuracy for each estimator, for a given context

s0v = s0;
idx_s0  = (s0v  == vidx);
% ==> to compute an accuracy value for ambiguous vertical orientation,
% choose a random 'cw' or 'ccw' groundtruth
s0v(idx_s0)   = c(binornd(1,0.5,sum(idx_s0),1)+1); 

% ==> accuracy computations
% ==> map accuracy
acc_map = (sum( (s0v < vidx)' & (map < vidx) ) + sum( (s0v > vidx)' & (map > vidx) )) / (sum(s0v < vidx) + sum(s0v > vidx));
% ==> mle accuracy
acc_mle = (sum( (s0v < vidx)' & (mle < vidx) ) + sum( (s0v > vidx)' & (mle > vidx) )) / (sum(s0v < vidx) + sum(s0v > vidx));    

% % ==> accuracy computations
% % ==> map accuracy
% acc_map = (sum( (s0v < vidx)' & (map < vidx) ) + sum( (s0v >= vidx)' & (map >= vidx) )) / (sum(s0v < vidx) + sum(s0v >= vidx));
% % ==> mle accuracy
% acc_mle = (sum( (s0v < vidx)' & (mle < vidx) ) + sum( (s0v >= vidx)' & (mle >= vidx) )) / (sum(s0v < vidx) + sum(s0v >= vidx));    

end