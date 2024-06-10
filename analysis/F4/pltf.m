function pltf(db, dp, aicd, llr, chpd, ppl, ctrs, ids)

% ==> marker colors
clr = [repmat([1,0.75,0],13,1); repmat([0.15,0.75,0.5],16,1)];
if (size(ppl,2) == 56) || (size(ppl,2) == 26) || (size(ppl,2) == 24) || (size(ppl,2) == 32) || (size(ppl,2) == 58)
    clr = [clr;clr];
end

msz = 60;

if strcmp(ctrs,'lo_ctr')
    alpha = repmat(0.25,29,1);
elseif strcmp(ctrs,'hi_ctr')
    alpha = ones(29,1);
elseif strcmp(ctrs,'lo_hi_ctr')
    alpha = [repmat(0.25,29,1); ones(29,1)];
end

% ==> choice predictivity and delta bias + delta uncertainty
[rb,pb] = corr(chpd(ppl), db(ppl));
[rp,pp] = corr(chpd(ppl), dp(ppl));

% ==> Use AIC computation from F2 instead of choice predictivity
[rb_aic,pb_aic] = corr(aicd(ppl), db(ppl));
[rp_aic,pp_aic] = corr(aicd(ppl), dp(ppl));

% ==> with LLR from F2
[rb_llr,pb_llr] = corr(llr(ppl), db(ppl));
[rp_llr,pp_llr] = corr(llr(ppl), dp(ppl));


% ==> pearson correlation between delta bias and delta uncertainty
% ==> ommitting session 1 here!
[rbp,pbp] = corr(db(ppl),dp(ppl));

% keyboard

subplot(3,3,1);
for i = ppl
    scatter(chpd(i), db(i), msz, 'o','filled','markeredgecolor', clr(i,:), 'markerfacecolor', clr(i,:), 'MarkerFaceAlpha', alpha(i),'MarkerEdgeAlpha', 1); 
    hold on; hold all;
    text(chpd(i)+0.005, db(i)-0.005, string([ids(i)]),'FontSize',4); 
end
title(['r = ',num2str(rb), ', p = ',num2str(pb)]);
xlabel('Choice Predictivity');
ylabel('\Delta bias');

subplot(3,3,2);
for i = ppl
    scatter(chpd(i), dp(i), msz, 'o','filled','markeredgecolor', clr(i,:), 'markerfacecolor', clr(i,:), 'MarkerFaceAlpha', alpha(i),'MarkerEdgeAlpha', 1); 
    hold on; hold all;
    text(chpd(i)+0.005, dp(i)-0.005,string([ids(i)]),'FontSize',4);
end
title(['r = ',num2str(rp), ', p = ',num2str(pp)]);
xlabel('Choice Predictivity');
ylabel('\Delta perceptual uncertainty (slope)');

subplot(3,3,3);
for i = ppl
    scatter(db(i), dp(i), msz, 'o','filled','markeredgecolor', clr(i,:), 'markerfacecolor', clr(i,:), 'MarkerFaceAlpha', alpha(i),'MarkerEdgeAlpha', 1); 
    hold on; hold all;
    text(db(i)+0.005,dp(i)-0.005,string([ids(i)]),'FontSize',4);
end
title(['r = ',num2str(rbp), ', p = ',num2str(pbp)]);
xlabel('\Delta bias')
ylabel('\Delta perceptual uncertainty (slope)')

subplot(3,3,4);
for i = ppl
    scatter(aicd(i), db(i), msz, 'o','filled','markeredgecolor', clr(i,:), 'markerfacecolor', clr(i,:), 'MarkerFaceAlpha', alpha(i),'MarkerEdgeAlpha', 1); 
    hold on; hold all;
    text(aicd(i)+0.005, db(i)-0.005, string([ids(i)]),'FontSize',4); 
end
title(['r = ',num2str(rb_aic), ', p = ',num2str(pb_aic)]);
xlabel('AIC \Delta');
ylabel('\Delta bias');

subplot(3,3,5);
for i = ppl
    scatter(aicd(i), dp(i), msz, 'o','filled','markeredgecolor', clr(i,:), 'markerfacecolor', clr(i,:), 'MarkerFaceAlpha', alpha(i),'MarkerEdgeAlpha', 1); 
    hold on; hold all;
    text(aicd(i)+0.005, dp(i)-0.005,string([ids(i)]),'FontSize',4);
end
title(['r = ',num2str(rp_aic), ', p = ',num2str(pp_aic)]);
xlabel('AIC \Delta');
ylabel('\Delta perceptual uncertainty (slope)');

subplot(3,3,7);
for i = ppl
    scatter(llr(i), db(i), msz, 'o','filled','markeredgecolor', clr(i,:), 'markerfacecolor', clr(i,:), 'MarkerFaceAlpha', alpha(i),'MarkerEdgeAlpha', 1); 
    hold on; hold all;
    text(llr(i)+0.005, db(i)-0.005, string([ids(i)]),'FontSize',4); 
end
title(['r = ',num2str(rb_llr), ', p = ',num2str(pb_llr)]);
xlabel('Log-likelihood ratio');
ylabel('\Delta bias');

subplot(3,3,8);
for i = ppl
    scatter(llr(i), dp(i), msz, 'o','filled','markeredgecolor', clr(i,:), 'markerfacecolor', clr(i,:), 'MarkerFaceAlpha', alpha(i),'MarkerEdgeAlpha', 1); 
    hold on; hold all;
    text(llr(i)+0.005, dp(i)-0.005,string([ids(i)]),'FontSize',4);
end
title(['r = ',num2str(rp_llr), ', p = ',num2str(pp_llr)]);
xlabel('Log-likelihood ratio');
ylabel('\Delta perceptual uncertainty (slope)');

end

