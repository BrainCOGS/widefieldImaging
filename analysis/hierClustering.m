function [clusID,dmat,lnk,CH] = hierClustering(weights,dmeasure,maxk,singKFlag)

% [clusID,dmat,lnk,CH] = hierClustering(weights,plotFlag,saveFlag,fnext)
%
% performs hiearchical clustering on matrix weights, of dimensions 
% observations (eg GLM weights) x nROI. The number of clusters is automatically 
% chosen by using the Calinski-Harabasz (CH) index to quantify "goodness of 
% clustering" (k = 1 to maxk) and determining the maximum CH. If a smaller k 
% has a CH value that is with 1 SD from the peak (calculated across k), it
% will be chosen instead. Linkage method is 'average'
% dmeasure is a string to choose distance measure, 'corr' for 1-r, where r
% is Pearson correlation, 'eucl' for eucledian distance
% singKFlag = 1 will assume maxk is optimal # clusters 
% 
% clusID is a nROI length vectors with cluster assignments, dmat is the 
% distance matrix used for clustering, lnk is linkage data for plotting dendograms 
% CH is the measure used to pick the number of clusters
%
% LP january 2014, modified april 2014, then april 2015

set(0,'RecursionLimit',1000)

% defaults / load variables
defaults={[];'eucl';10;0};
inputnames = {'weights';'dmeasure';'maxk';'singKFlag'};
if nargin < length(defaults)
    for i = nargin+1:length(defaults)
        temp=defaults{i};
        eval([inputnames{i} '=temp;']);
    end
end

if isempty(weights)
    error('I need at least one input')
end

switch dmeasure
    case 'correlation'
        rDist = pdist(weights,'correlation'); % 1-r, where r is Pearson correlation
        dmat = squareform(rDist);
        lnk = linkage(rDist,'average'); % linkage matrix
    case 'eucl'
        rDist = pdist(weights); % use euclidean distance 
        dmat = squareform(rDist);
        lnk = linkage(rDist,'ward'); % linkage matrix
end
% linkage stats
inc = inconsistent(lnk);
cp = cophenet(lnk,rDist);

% % plot dendrogram 
% figure;
% dendrogram(lnk); xlabel('ROI #'); ylabel('distance')
% % close
% % saveas(gcf,['hc_dendrogram' fnext]); 
% pause; close

% run 1 to maxk clusters and quantify goodness of clustering to decide
% optimal k. uses the Calinski-Harabasz index.   

if ~singKFlag
    cid = zeros(maxk,size(weights,1)); ks = 1:maxk;
    for i = ks
        cid(i,:) = cluster(lnk,'maxclust',i);
    end
    CH = hcval(weights,cid');
    maxclusters = ks(find(CH==max(CH),1,'first'));
%     
%     maxclusters = min(find(CH>=(max(CH)-std(CH(2:end)))));
%     
%     % automatically pick first peak CH value that is above maximum - std of CH
%     % across different k
%     figure; hold on
%     plot(CH)
%     plot([0 maxk],[max(CH)-std(CH(2:end)) max(CH)-std(CH(2:end))],'r--')
%     plot(maxclusters,CH(maxclusters),'ro')
%     maxclusters = input('number of clusters: '); % manually pick first peak CH value
%     close
    
    % final clustering with optimal k
    clusID = cid(maxclusters,:)';
else
    clusID = cluster(lnk,'maxclust',maxk);
    CH     = [];
    
end

% try
%     maxclusters = find(CH>=std(CH(1:end)),1,'first');
% catch
%     maxclusters = maxk;
% end
%saveas(gcf,['CH' fnext]); 


