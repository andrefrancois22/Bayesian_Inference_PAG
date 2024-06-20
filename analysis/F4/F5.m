function F5(db,dp)
    % INPUT: 
    % (1) db, dp (delta bias and delta uncertainty)
    
    % ==> helper functions
    % ==> fisher transformation function
    z = @(r) (1/2)*log((1 + r)/(1 - r));
    % ==> inverse of the z transformation
    zi = @(z) (exp(2*z) - 1)/(exp(2*z) + 1);
    % ==> normalization function
    normf = @(x) (x - mean(x))/std(x);

    % => Monkey F
    rFl = corr(db(1:13,1),dp(1:13,1));
    rFh = corr(db(1:13,2),dp(1:13,2));
    % => Monkey J
    rJl = corr(db(14:end,1),dp(14:end,1));
    rJh = corr(db(14:end,2),dp(14:end,2));    

    % ==> error bounds for plots in panel
    % => standard deviation (1/sqrt(N - 3))
    % => monkey F (N = 13), so SD denom = sqrt(10)
    sdFl_ub = z(rFl) + 1/sqrt(10); sdFl_lb = z(rFl) - 1/sqrt(10); 
    sdFh_ub = z(rFh) + 1/sqrt(10); sdFh_lb = z(rFh) - 1/sqrt(10);
    % => monkey J (N = 16), so SD denom = sqrt(13)
    sdJl_ub = z(rJl) + 1/sqrt(13); sdJl_lb = z(rJl) - 1/sqrt(13); 
    sdJh_ub = z(rJh) + 1/sqrt(13); sdJh_lb = z(rJh) - 1/sqrt(13);

    % ~~~~~~~~~ NOTE: (errorbar(X,Y,NEG(i),POS(i))) ~~~~~~~~~
    % ==> errorbar native plot function allows NEG(i), POS(i), which are ~distances~ above and below
    % correlations for drawing error bars in the plot. So these have to be computed as follows:

    % => NEG(i) is correlation r - z^(-1)(SD_LB), where SD_LB is the lower
    % bound correlation in z coordinates, e.g. 
    % ==> SD_LB = z(r) - 1/sqrt(N-3)
    % so the ~distance~ below r for NEG(i) is r - z^(-1)(SD_LB)
    % (actual correlation minus the lower bound correlation)

    % => POS(i) is correlation z(-1)(SD_UB), where SD_UB is the higher 
    % bound correlation in z coordinates, e.g.
    % ==> SD_UB = r(z) + 1/sqrt(N-3)
    % so the ~distance~ above r for POS(i) is z(-1)(SD_UB) - r
    % (upper bound correlation minus actual correlation)

    % => Monkey colors
    mclrs = {[1,0.75,0],[0.15,0.75,0.5]};    
    figure; set(gcf,'color','white');
    hold on; hold all;
    % ==> monkey F 
    errorbar([4,5], ...
             [rFl, rFh], ... % ==> original correlations
             [rFl-zi(sdFl_lb), rFh-zi(sdFh_lb)],[zi(sdFl_ub)-rFl, zi(sdFh_ub)-rFh] , ...
             's','color',mclrs{1});
    % ==> monkey J
    errorbar([2,3], ...
             [rJl, rJh], ... % ==> original correlations
             [rJl-zi(sdJl_lb), rJh-zi(sdJh_lb)], [zi(sdJl_ub)-rJl, zi(sdJh_ub)-rJh] , ...
             's','color',mclrs{2}); 
    xlim([0,7]);
    ylim([-0.1,1]);

    % ==> z-scores across animals
    dbz = nan(29,2);
    dpz = nan(29,2);

    % ==> low & high contrast
    for c = 1:2
        % ==> animal F
        dbz(1:13,c) = normf(db(1:13,c));
        dpz(1:13,c) = normf(dp(1:13,c));
        % ==> animal J
        dbz(14:29,c) = normf(db(14:29,c));
        dpz(14:29,c) = normf(dp(14:29,c));
    end

    % ==> correlation
    [r_z_lcr, ~] = corr(dbz(:,1),dpz(:,1));
    [r_z_hcr, ~] = corr(dbz(:,2),dpz(:,2));

    % ==> correlation across all contrasts and animals
    [rz_dp, ~] = corr([dbz(:,1); dbz(:,2)],[dpz(:,1); dpz(:,2)]);

    % ==> lower and upper bound values for SE
    % => low contrast cases
    lcr_se_lb = z(r_z_lcr) - 1/sqrt(26); 
    lcr_se_ub = z(r_z_lcr) + 1/sqrt(26);
    % => high contrast cases
    hcr_se_lb = z(r_z_hcr) - 1/sqrt(26); 
    hcr_se_ub = z(r_z_hcr) + 1/sqrt(26);
    % => all data
    se_lb = z(rz_dp) - 1/sqrt(29*2 - 3); 
    se_ub = z(rz_dp) + 1/sqrt(29*2 - 3);
    % ==> draw plot panel
    errorbar([6,7,1], ...
             [r_z_lcr, r_z_hcr, rz_dp], ...
             [r_z_lcr-zi(lcr_se_lb), r_z_hcr-zi(hcr_se_lb), rz_dp-zi(se_lb)], ... 
             [zi(lcr_se_ub)-r_z_lcr, zi(hcr_se_ub)-r_z_hcr, zi(se_ub)-rz_dp], ...
             'ks');

    % => axis labels
    ax = gca;      
    ax.XTickLabel = {[],'all','J-lo','J-hi','F-lo','F-hi','lo','hi'}; 
    xlabel('Monkey and Stimulus Contrast Condition')  
    ylabel('Association (r)');    
             
end

