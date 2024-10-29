function [acc_map, acc_mle] = calcaccf(s0,vidx,map,mle,c)

% ~~~~~~~~~~~~~~~ ACCURACY CALCULATION FUNCTION ~~~~~~~~~~~~~~~
% This calculates the accuracy for each estimator, for a given context
% ~~~~~~~~~~~~ INPUT ~~~~~~~~~~~~
% s0             ==> physical stimulus value for each simulated trial
% vidx           ==> index for simulated neutral (vertical) stimulus
% map            ==> maximum a posteriori estimates - orientation values 
% mle            ==> maximum likelihood estimates - orientation values
% c              ==> [1,7]: a ccw and a cw value (for handling ambiguous cases with 7 total orienations)
% ~~~~~~~~~~~~ OUTPUT ~~~~~~~~~~~~
% acc_mle        ==> MLE accuracy
% acc_map        ==> MAP accuracy


s0v = s0;
idx_s0  = (s0v  == vidx);
% ==> to compute an accuracy value for ambiguous vertical orientation,
% choose a random 'cw' or 'ccw' groundtruth
s0v(idx_s0)   = c(binornd(1,0.5,sum(idx_s0),1)+1); 

% % ==> accuracy computations
% % ==> map accuracy
% acc_map = (sum( (s0v < vidx)' & (map < vidx) ) + sum( (s0v > vidx)' & (map > vidx) )) / (sum(s0v < vidx) + sum(s0v > vidx));
% % ==> mle accuracy
% acc_mle = (sum( (s0v < vidx)' & (mle < vidx) ) + sum( (s0v > vidx)' & (mle > vidx) )) / (sum(s0v < vidx) + sum(s0v > vidx));    

% ==> accuracy computations
% ==> map accuracy
acc_map = (sum( (s0v < vidx)' & (map < vidx) ) + sum( (s0v >= vidx)' & (map >= vidx) )) / (sum(s0v < vidx) + sum(s0v >= vidx));
% ==> mle accuracy
acc_mle = (sum( (s0v < vidx)' & (mle < vidx) ) + sum( (s0v >= vidx)' & (mle >= vidx) )) / (sum(s0v < vidx) + sum(s0v >= vidx));    


end
