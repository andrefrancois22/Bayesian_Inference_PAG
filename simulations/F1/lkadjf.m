function lk = lkadjf(p,lk,sig,K,s,d,vidx,pd)
% this function sums truncated density values in likelihood 

% ==> padded orientation range (to capture full density tails)
d2 = -pd:pd;
lk_tr = p(s',sig,d2);

% ==> range must be large enough for normal pdf to decay fully (avoid truncation)
assert( (sum(sum(lk_tr(:,[1,end]))) == 0), '>>>>>> padding parameter "pd" is too small! <<<<<<' )

% ==> strfind works here
% => get indices of orientations in d inside larger d2 range
idx = strfind(d2, d):(strfind(d2, d)+length(d)-1);

% ==> create mask for sum
m = ones(length(d2),1);
m(idx)=0;
% => mask array
msk = repmat(m',[K,1]);
% => get truncated density values
dn = lk_tr.*msk; 

% ==> sum truncated values on each side of vertical orientation
% => this yields a K x 2 matrix containing sum of truncated values on either end
dns = [sum(dn(:,1:(idx(vidx)-1)),2),sum(dn(:,(idx(vidx)+1):end),2)];

% ==> add truncated density to edges of likelihood functions for all k
lk(:,1)   = lk(:,1)   + dns(:,1);
lk(:,end) = lk(:,end) + dns(:,2);
% ==> normalize result
lk = lk./sum(lk,2);

end

