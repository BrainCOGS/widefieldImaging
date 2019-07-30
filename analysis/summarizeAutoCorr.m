function acorrSumm = summarizeAutoCorr(mice)

% decodeSumm = summarizeDecoding(mice,cfg,whichDecoder)
% summarizes pairwise decoding weight correlations as a function of
% distance
%         mice: cell array with mouse name list
%          cfg: structure with analysis config, pass empty to load defaults    
% whichDecoder: 'evidence', 'choice', 'prevchoice' or 'viewangle'

%% initialize
if nargin < 1 || isempty(mice); mice           = {'ai2';'ai3';'ai5';'ai7';'ai9';'ai10'}; end

try

%% compile from saved data files
tic
fprintf('compiling data...\n')
recls = cell(numel(mice),1);
for iMouse = 1:numel(mice)
  
  fprintf('\tmouse %d / %d',iMouse,numel(mice))
  
  [recls{iMouse},rootdir] = getFullRecPath(mice{iMouse},isThisSpock);
  for iRec = 1:numel(recls{iMouse})
    fprintf('.')
    cd(recls{iMouse}{iRec})
    try
      load('pxlAutoCorr.mat','dffAutoCorr'); 
      acorrSumm.mouse(iMouse).corrByDist(:,iRec) = dffAutoCorr.corrByDist_mean;
      acorrSumm.distAxis                         = dffAutoCorr.distAxis;
      acorrSumm.mouse(iMouse).acorr(:,iRec)      = dffAutoCorr.acorrs_mean;
      acorrSumm.timeAxis_acorr                   = dffAutoCorr.taxis;
    catch
      acorrSumm.mouse(iMouse).corrByDist(:,iRec) = nan(numel(acorrSumm.distAxis),1);
      acorrSumm.mouse(iMouse).acorr(:,iRec)      = nan(numel(acorrSumm.timeAxis_acorr),1);
    end
  end
  
  %% mouse average
  acorrSumm.corrByDist(:,iMouse) = nanmean(acorrSumm.mouse(iMouse).corrByDist,2);
  acorrSumm.acorr(:,iMouse)      = nanmean(acorrSumm.mouse(iMouse).acorr,2);
  fprintf('\n')
  
end

%% save
cd(rootdir)
fn = 'autoCorrSummary'; 
save(fn,'acorrSumm')

%% plot
wf          = widefieldParams;
figure; wf.applyFigDefaults(gcf,[2 1],'w');

subplot(1,2,1); hold on
xaxis = acorrSumm.timeAxis_acorr;
thism = mean(acorrSumm.acorr,2);
thiss = std(acorrSumm.acorr,0,2)./sqrt(numel(mice)-1);
errorshade(thism,thiss,[0 0 0],[.7 .7 .7],xaxis,1.5);
yl = get(gca,'ylim');
plot([0 0],yl,'--','color',[.5 .5 .5]);
wf.applyAxisDefaults(gca,'k');
wf.applyAxisLbls(gca,'lag (s)','auto-corr','Pxl autocorr')

subplot(1,2,2); hold on
xaxis = acorrSumm.distAxis;
thism = mean(acorrSumm.corrByDist,2);
thiss = std(acorrSumm.corrByDist,0,2)./sqrt(numel(mice)-1);
errorshade(thism,thiss,[0 0 0],[.7 .7 .7],xaxis,1.5);
wf.applyAxisDefaults(gca,'k');
wf.applyAxisLbls(gca,'Dist (mm)','corr','Spatial autocorr')

saveas(gcf,[fn '.png'],'png')

catch ME
  displayException(ME)
end
end