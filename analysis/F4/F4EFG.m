function [pJl,pFl,pJh,pFh,pl,ph,pall] = F4EFG(d, aic_d, metric)
    % INPUT: 
    % (1) d (delta metric data, can be delta bias, or delta uncertainty)
    % ==> delta bias is the difference betweem the pairs of decision
    % criteria estimated for the PFs from the median split DV dynamic range analysis 
    % (2) aic_d (delta AIC GLM model comparison, computed for F2)
    % (3) metric (just a string to indicate which metric was passed)
    % for db (delta bias), string should be 'bias'
    % for dp (delta uncertainty), string should be 'uncertainty'

    % OUTPUT:
    % p-values for wilcoxon tests
    
    % ==> helper functions (for association r value panel)
    % ==> fisher transformation function
    z = @(r) (1/2)*log((1 + r)/(1 - r));
    % ==> inverse of the z transformation
    zi = @(z) (exp(2*z) - 1)/(exp(2*z) + 1);

    % ==> plot average delta Bias
    figure; set(gcf,'color','white'); set(gcf,'Position',[675 560 1014 402]);
    subplot(1,2,1);
    hold on; hold all;
    % ==> monkey F - lo & hi contrast
    errorbar([1,2], ...
             [mean(d(1:13,1)), mean(d(1:13,2))], ...
             [std(d(1:13,1))/sqrt(13), std(d(1:13,2))/sqrt(13)] , ...
             'o','color','k');
    % ==> monkey J - lo & hi contrast
    errorbar([3,4], ...
             [mean(d(14:end,1)), mean(d(14:end,2))], ...
             [std(d(14:end,1))/sqrt(16), std(d(14:end,2))/sqrt(16)] , ...
             'ko');     
    % ==> lo & hi contrast     
    errorbar([5,6], ...
             [mean(d(:,1)), mean(d(:,2))], ...
             [std(d(:,1))/sqrt(length(d(:,1))), std(d(:,2))/sqrt(length(d(:,2)))] , ...
             'ko');
    % ==> all
    errorbar([7], ...
             [mean([d(:)])], ...
             [std([d(:)])/sqrt(length(d(:))*2)], ...
             'ko');
    xlim([0,7]);
    % => axis labels
    ax = gca;      
    title(['\Delta ', metric, ' results'])
    ylabel(['\Delta ', metric])
    ax.XTickLabel = {[],'F-lo','F-hi','J-lo','J-hi','lo','hi','all'}; 
    xlabel('Monkey and Stimulus Contrast Condition')    
    
    % ==> wilcoxon tests
    % => low contrast cases
    [pJl,~] = signrank(d(14:end,1));
    [pFl,~] = signrank(d(1:13,1));
    % => high contrast cases
    [pJh,~] = signrank(d(14:end,2));
    [pFh,~] = signrank(d(1:13,2));
    % => high & low separately
    [pl,~] = signrank(d(:,1));
    [ph,~] = signrank(d(:,2));
    % => all
    [pall,~] = signrank([d(:)]);

    % => Monkey F
    rFl = corr(d(1:13,1),aic_d(1:13));
    rFh = corr(d(1:13,2),aic_d(1:13));
    % => Monkey J
    rJl = corr(d(14:end,1),aic_d(14:end));
    rJh = corr(d(14:end,2),aic_d(14:end));

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

    subplot(1,2,2);
    hold on; hold all;
    % ==> monkey F 
    errorbar([1,2], ...
             [rFl, rFh], ... % ==> original correlations
             [rFl-zi(sdFl_lb), rFh-zi(sdFh_lb)],[zi(sdFl_ub)-rFl, zi(sdFh_ub)-rFh] , ...
             'o','color','k');
    % ==> monkey J
    errorbar([3,4], ...
             [rJl, rJh], ... % ==> original correlations
             [rJl-zi(sdJl_lb), rJh-zi(sdJh_lb)], [zi(sdJl_ub)-rJl, zi(sdJh_ub)-rJh] , ...
             'ko'); 

    % ==> z-scores across animals
    dz = nan(29,2);

    % ==> Animal F
    % ==> low contrast
    dz(1:13,1) = (d(1:13,1) - mean(d(1:13,1)))/std(d(1:13,1));
    % ==> high contrast
    dz(1:13,2) = (d(1:13,2) - mean(d(1:13,2)))/std(d(1:13,2)); 
    % ==> Animal J
    % ==> low contrast
    dz(14:29,1) = (d(14:29,1) - mean(d(14:29,1)))/std(d(14:29,1));
    % ==> high contrast
    dz(14:29,2) = (d(14:29,2) - mean(d(14:29,2)))/std(d(14:29,2)); 

    % ==> standardize aic by monkey
    aic_d_z(1:13)  = (aic_d(1:13)  - mean(aic_d(1:13)))/std(aic_d(1:13));
    aic_d_z(14:29) = (aic_d(14:29) - mean(aic_d(14:29)))/std(aic_d(14:29));
    
    % ==> correlation
    [r_z_lcr, ~] = corr(aic_d_z',dz(:,1));
    [r_z_hcr, ~] = corr(aic_d_z',dz(:,2));

    % ==> correlation across all contrasts and animals
    [rz_d, ~] = corr([aic_d_z,aic_d_z]',[dz(:,1); dz(:,2)]);

    % ==> Estimate standard errors with z transformation and inverse z^(-1)
    % => plot SE lower bound (lb) and higher bound (hb) for low contrast cases
    lcr_se_lb = z(r_z_lcr) - 1/sqrt(26); 
    lcr_se_ub = z(r_z_lcr) + 1/sqrt(26);
    % => plot SE lower bound (lb) and higher bound (hb) for high contrast cases
    hcr_se_lb = z(r_z_hcr) - 1/sqrt(26); 
    hcr_se_ub = z(r_z_hcr) + 1/sqrt(26);
    % => all
    se_lb = z(rz_d) - 1/sqrt(29*2 - 3); 
    se_ub = z(rz_d) + 1/sqrt(29*2 - 3);

    % => draw panel
    errorbar([5,6,7], ...
             [r_z_lcr, r_z_hcr, rz_d], ...
             [r_z_lcr-zi(lcr_se_lb), r_z_hcr-zi(hcr_se_lb), rz_d-zi(se_lb)], ...
             [zi(lcr_se_ub)-r_z_lcr, zi(hcr_se_ub)-r_z_hcr, zi(se_ub)-rz_d], ...
             'ro');
    xlim([0,7]);
    ylim([0,1]); 
    
    % => axis labels
    ax = gca;    
    title(['\Delta ', metric, ' results'])    
    ylabel(['\Delta ', metric, ' and \Delta AIC - Association (r) value'])
    ax.XTickLabel = {[],'F-lo','F-hi','J-lo','J-hi','lo','hi','all'};   
    xlabel('Monkey and Stimulus Contrast Condition')

end

