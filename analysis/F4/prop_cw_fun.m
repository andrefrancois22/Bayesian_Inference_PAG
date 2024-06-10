
% ==> low dynamic range cases

% ==> compute ratios for actual data
num_cw = sum(cw_lcr_cng & ori==or(th) & ldyn); 
denom_cw = num_cw + sum(cw_lcr_icg & ori==or(th) & ldyn);
% ==> proportion cw choice
prop_cw_lcr_ldyn(th) = num_cw/denom_cw;

num_ccw = sum(ccw_lcr_icg & ori==or(th) & ldyn);
denom_ccw = num_ccw + sum(ccw_lcr_cng & ori==or(th) & ldyn);
% ==> proportion cw choice
prop_ccw_lcr_ldyn(th) = num_ccw/denom_ccw;   

% ==> compute ratios for actual data
num_cw = sum(cw_hcr_cng & ori==or(th) & ldyn);
denom_cw = num_cw + sum(cw_hcr_icg & ori==or(th) & ldyn);
% ==> proportion cw choice
prop_cw_hcr_ldyn(th) = num_cw/denom_cw;

num_ccw = sum(ccw_hcr_icg & ori==or(th) & ldyn);
denom_ccw = num_ccw + sum(ccw_hcr_cng & ori==or(th) & ldyn);
% ==> proportion cw choice
prop_ccw_hcr_ldyn(th) = num_ccw/denom_ccw; 

% ==> high dynamic range cases

% ==> compute ratios for actual data
num_cw = sum(cw_lcr_cng & ori==or(th) & hdyn); 
denom_cw = num_cw + sum(cw_lcr_icg & ori==or(th) & hdyn);
% ==> proportion cw choice
prop_cw_lcr_hdyn(th) = num_cw/denom_cw;

num_ccw = sum(ccw_lcr_icg & ori==or(th) & hdyn);
denom_ccw = num_ccw + sum(ccw_lcr_cng & ori==or(th) & hdyn);
% ==> proportion cw choice
prop_ccw_lcr_hdyn(th) = num_ccw/denom_ccw;   

% ==> compute ratios for actual data
num_cw = sum(cw_hcr_cng & ori==or(th) & hdyn);
denom_cw = num_cw + sum(cw_hcr_icg & ori==or(th) & hdyn);
% ==> proportion cw choice
prop_cw_hcr_hdyn(th) = num_cw/denom_cw;

num_ccw = sum(ccw_hcr_icg & ori==or(th) & hdyn);
denom_ccw = num_ccw + sum(ccw_hcr_cng & ori==or(th) & hdyn);
% ==> proportion cw choice
prop_ccw_hcr_hdyn(th) = num_ccw/denom_ccw; 
