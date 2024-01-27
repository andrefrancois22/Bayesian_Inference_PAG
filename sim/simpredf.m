function [mle,map,pos,acc_mle,acc_map,lk,po] = simpredf(pr,s0,vidx,sig,n_o,K,LK_TRC_FLAG)
    % Computes likelihood, posterior, MAP, pdf density >= cw orientation and MLE accuracy

    % => delta away from stimulus s can be one of the following:
    % { -3, -2, -1, 0, +1, +2, +3 }
    d = (1:n_o) - vidx;

    % => define standard normal pdf function
    p = @(mu,sig,ds) normpdf(ds,mu-vidx,sig);

    % => get probabilities of deviations centered on each kth ground truth stimulus (from s0)
    p_eps = p(s0',sig,d);
    % => normalization necessary
    p_eps = p_eps./sum(p_eps,2); 

    % ==> (3) perceptual noise step
    % ==> random epsilon displacement angle index from ground truth orientation for each kth stimulus s
    % s are the sensory stimulus orientations (s0 ground truth perceived with some added sensory noise)
    s = arrayfun(@(i) randsample(1:n_o,1,true,p_eps(i,:)), 1:K);

    % ==> (4) define likelihood
    lk = p(s',sig,d); %**** note edge cases - truncated values
    % => normalize 
    lk = lk./sum(lk,2); 
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ==> special option: sum truncated density values and add to extremes
    if ~LK_TRC_FLAG
        pd = 1000;
        lk = lkadjf(p,lk,sig,K,s,d,vidx,pd);
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % ==> (5) inference step - posterior
    po = repmat(pr,[K,1]).*lk ./ dot(repmat(pr,[K,1])',lk')';    

    % ==> compute MLE, MAP & mass orientation estimates  
    % [~,mle] = max(lk,[],2); %**** unnecessary - just s
    mle = s';
    
    % ==> map estimate
    %[~,map] = max(po,[],2); 
    % ==> posterior log likelihood ratio
    map = nan(K,1);
    indChooseCW  = log(sum(po(:,(vidx+1):n_o), 2)   ./sum(po(:,1:(vidx-1)), 2)) > 0;
    indChooseCCW = log(sum(po(:,(vidx+1):n_o), 2)   ./sum(po(:,1:(vidx-1)), 2)) < 0;
    indGuess     = log(sum(po(:,(vidx+1):n_o), 2)   ./sum(po(:,1:(vidx-1)), 2)) == 0;
    map(indChooseCW) = n_o; 
    map(indChooseCCW) = 1;
    map(indGuess) = datasample([1 n_o], sum(indGuess), 'replace', true);   

    % ==> posterior sample
    pos = arrayfun(@(i) randsample(1:n_o,1,true,po(i,:)), 1:K)';

    %**** if a response in mle or map is a (vertical) '4' ~> then make a random 'cw' or 'ccw' choice         
    % ~~~~~~~~~~~~~~ if response in mle or map is a vertical orientation ~~~~~~~~~~~~~~~
    % ==> find indices of 'vertical' response
    idx_mle = (mle == vidx);
    idx_map = (map == vidx);
    c=[1,n_o];
    % ==> choose randomly between 'cw' and 'ccw' orientation ('1' or '7' coin toss)    
    mle(idx_mle) = c(binornd(1,0.5,sum(idx_mle),1)+1); 
    map(idx_map) = c(binornd(1,0.5,sum(idx_map),1)+1); 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

    % ==> checks
    assert(isempty(mle(mle == vidx)));
    assert(isempty(map(map == vidx)));    
    
    % ==> calculate accuracy for each estimator and given a context
    [acc_map, acc_mle] = calcaccf(s0,vidx,map,mle,c);       
    
    % ==> terminal updates...
    fprintf('MAP acc=%d, MLE acc=%d for K=%d sample trials, sigma=%d...\n',acc_map,acc_mle,K,sig)

% keyboard

end

