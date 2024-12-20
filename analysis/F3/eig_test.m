% ==> draw figure
fg = figure(); fg.Position = [328 68 970 892]; clf; set(gcf,'Color','w');
% ==> subplot indices
sbpi = [1,3;2,4];

% ==> bootstrap sample number
BT = 1000;
    
for m = 1:length(M)

    % ==> early & late window dv values 
    v = vertcat(dvw{rgs{m},:,:});
    % ==> for random sample with replacement
    N = length(v);
    K = N;

    % ==> store slopes
    s = nan(1,BT);

    % ==> sample
    for b = 1:BT

        % ==> indices for random sample with replacement
        idx = randsample(N,K,true);
        % ==> covariance matrix
        cvb = cov(v(idx,:));
        % ==> eigendecomposition
        [V,D] = eig(cvb);
        % ==> get index of max eigenvalue
        [~,I] = max(D);
        % ==> mean of the distribution
        mub = mean(cvb)';

        % ==> angle of first eigenvector
        [TH,~] = cart2pol(V(1,I(2)),V(2,I(2)));
        % ==> store slope
        s(b) = TH;
    end

    % ==> get confidence intervals - is 45-deg inside that interval for s?
    % or 0 for s - pi/4
    prc = prctile((s - pi/4),[2.5 97.5]);
    fprintf('Percentile values for confidence interval for 2.5 and 97.5 density are %d and %d...\n',prc(1),prc(2))
    % ==> check that 0 is not within the confidence interval - e.g. prc
    % values are larger than 0
    assert(unique(prc > 0))
    if unique(prc > 0)
        fprintf('Confidence interval bounds for (s - pi/4) are larger than 0...\n')
    else
        fprintpf('Confidence interval bounds for (s - pi/4) include 0...\n')
    end

    % ==> compute difference between estimated slopes and 'null' slope (e.g. x=y, or TH = pi/4)
    ds = s - pi/4;
    % ==> test that ds deviated significantly from zero (means slope angle is greater than pi/4 e.g. x=y)
    %[ps, ~] = signrank(ds);
    % ==> no good - signrank inflates significance because of N

    fprintf('first eigenvector slope for monkey %s...\n',M{m});
    
    % ==> subplot
    subplot(2,2,sbpi(1,m));
    scatter(v(:,1),v(:,2),'ko','filled');
    hold on; hold all;
    for xc = 1:length(cx)
        for rc = 1:length(cr)
            % ==> scatter early and late dv values
            v = vertcat(dvw{rgs{m},xc,rc});
            scatter(v(:,1),v(:,2), 140, 'o','markerfacecolor', clrs{m}, 'markeredgecolor', 'w');
        end
    end
    % ==> lines
    plot(-0.6:0.1:0.6,-0.6:0.1:0.6,'r--')
    plot(-0.6:0.1:0.6,zeros([1,length(-0.6:0.1:0.6)]),'k--')
    plot(zeros([1,length(-0.6:0.1:0.6)]),-0.6:0.1:0.6,'k--')
    axis equal; axis square;
    xlabel('Average (-800ms to -600ms) window');
    ylabel('Average (-500ms to -300ms) window');
    title(['First PC slope (deg): ', num2str(rad2deg(THs(m)))])
    % ==> original eigenvectors
    plot([MUs(1,m),Vs(1,m)],[MUs(2,m),Vs(2,m)],'r-','linewidth', 3)   
    xlim(xyl{m}); 
    ylim(xyl{m});   
    drawnow;

    % ==> plot distribution of PC angles
    bns = 25;
    subplot(2,2,sbpi(2,m)); cla;
    h = histogram(s,bns,'Facecolor','r'); hold on; hold all;
    xlabel('First PC slope (radians)');
    ylabel('frequency');
    title(['distribution of slopes of first PC: ',M{m}]);
    plot(repmat(pi/4,[1,length(0:max(h.Values))]),0:max(h.Values),'r--')
    text(pi/4 + 0.025,10,'\pi/4 equality line (x=y) angle (slope)','Rotation',90,'Color','r')
    xlim([pi/5,pi/2]);

end