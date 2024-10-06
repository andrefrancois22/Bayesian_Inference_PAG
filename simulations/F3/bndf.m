function [csens_cx_bnd] = bndf(csens_cx,bnd,N, M_FLAG)

% ==> sets DV values that reach a bound to the bound value, including all
% subsequent DV values following the bound crossing

if strcmp(M_FLAG, 'BOUND')   
    % ==> set DV values for timepoints over trials that reach or exceed a bound
    % to the bound value.
    csens_cx(csens_cx >=  bnd) =  bnd;
    csens_cx(csens_cx <= -bnd) = -bnd;
    
    % ==> set all dv values beyond the bound to the bound value
    % ==> positve bound
    for i = 1:N
        bidx = find(csens_cx(i,:) >=  bnd, 1, 'first');
        % => set to bound
        if ~isempty(bidx)
            csens_cx(i,bidx:end) = bnd;
        end
    end
    % ==> negative bound
    for i = 1:N
        bidx = find(csens_cx(i,:) <=  -bnd, 1, 'first');
        % => set to bound
        if ~isempty(bidx)
            csens_cx(i,bidx:end) = -bnd;
        end
    end
end

% ==> output
csens_cx_bnd = csens_cx;



