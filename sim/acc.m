function [mle,map,pos,acc_mle,acc_map,acc_dns,acc_pos,lk,po] = simpredf(p_cw,s0,vidx,sig,n_o,K)
    % Computes likelihood, posterior, MAP, pdf density >= cw orientation and MLE accuracy

    % => delta away from stimulus s can be one of the following:
    % { -3, -2, -1, 0, +1, +2, +3 }
    d = (1:n_o) - vidx;

    % => define standard normal pdf function
    p = @(mu,sig) normpdf(d,mu-vidx,sig);

    % => get probabilities of deviations centered on each kth ground truth stimulus (from s0)
    p_eps = p(s0',sig);
    % => normalization appears necessary
    p_eps = p_eps./sum(p_eps,2); 

    % ==> (3) perceptual noise step
    % ==> random epsilon displacement angle index from ground truth orientation for each kth stimulus s
    % s are the sensory stimulus orientations (s0 ground truth perceived with some added sensory noise)
    s = arrayfun(@(i) randsample(1:n_o,1,true,p_eps(i,:)), 1:K);

    % ==> (4) define likelihood
    lk = p(s',sig); %**** note edge cases - truncated values
    % => normalize 
    lk = lk./sum(lk,2); 

    % ==> (5) inference step - posterior
    po = repmat(p_cw,[K,1]).*lk ./ dot(repmat(p_cw,[K,1])',lk')';

    % ==> compute MLE, MAP & mass orientation estimates  
    mle = s';
    % [~,mle] = max(lk,[],2); %**** unnecessary - just s
    % assert(sum(mle==s') == K)
    % ==> map estimate
    [~,map] = max(po,[],2);
    % ==> posterior sample
    pos = arrayfun(@(i) randsample(1:n_o,1,true,po(i,:)), 1:K)';

    %**** evaluate 0 (vertical case randomly
    %**** if a response in mle or map is a '4' ~> random 'cw' or 'ccw' choice         
    % ~~~~~~~~~~~~~~ if response in mle or map is a vertical orientation ~~~~~~~~~~~~~~~
    % ==> find indices of 'vertical' response
    idx_mle = (mle == vidx);
    idx_map = (map == vidx);
    % ==> choose randomly between 'cw' and 'ccw' orientation ('1' or '7')
    mle(idx_mle) = randi([1,n_o],[sum(idx_mle),1])';
    map(idx_map) = randi([1,n_o],[sum(idx_map),1])';
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    
%     % ==> mle accuracy
%     acc_mle  = sum( (s0 > vidx)' & (mle > vidx) ) / sum(s0 > vidx); 
%     % ==> map accuracy
%     acc_map  = sum( (s0 > vidx)' & (map > vidx) ) / sum(s0 > vidx);
%     % ==> density
%     acc_dns  = sum( (s0 > vidx)' & (sum(po(:,(vidx+1):end),2) > 0.5)) / sum(s0 > vidx);
%     % ==> posterior sample
%     acc_pos  = sum( (s0 > vidx)' & (pos > vidx) ) / sum(s0 > vidx);
%        
    % ==> mle accuracy
    acc_mle  = sum( (s0 < vidx)' & (mle < vidx) ) / sum(s0 < vidx); 
    % ==> map accuracy
    acc_map  = sum( (s0 < vidx)' & (map < vidx) ) / sum(s0 < vidx);
    % ==> density
    acc_dns  = sum( (s0 < vidx)' & (sum(po(:,1:vidx),2) > 0.5)) / sum(s0 < vidx);
    % ==> posterior sample
    acc_pos  = sum( (s0 < vidx)' & (pos < vidx) ) / sum(s0 < vidx);
           
    
    % ==> terminal updates...
    fprintf('MAP acc=%d, MLE acc=%d, mass acc=%d for K=%d sample trials, sigma=%d...\n',acc_map,acc_mle,acc_dns,K,sig)
end

