function seq = activitySequences(avgDff,cfg,plotFlag)
  
% seq = activitySequences(avgDff,cfg)
% compares pxl(roi) activity sequences by choice and maze 
% avgDff is the output structure of avgDffByTrialEpoch() or avgDffROI()
% cfg is optional analysis configuration structure
% plotFlag is boolean to plot analysis results

if nargin < 2; cfg.maxBinSequence = 250;  end
if nargin < 3; plotFlag           = true; end

%% get left and right choice avgs for each maze
Raccumul  = avgDff.maze(end).correctR.space.frames;
Laccumul  = avgDff.maze(end).correctL.space.frames;
RvisGuide = avgDff.maze(1).correctR.space.frames;
LvisGuide = avgDff.maze(1).correctL.space.frames;

%% reshape if in pxl format
if size(Raccumul,3) > 1
  Raccumul  = conditionDffMat(Raccumul);
  Laccumul  = conditionDffMat(Laccumul);
  RvisGuide = conditionDffMat(RvisGuide);
  LvisGuide = conditionDffMat(LvisGuide);
end

%% exclude later bins if desired
bins      = avgDff.spaceBins;
isGoodBin = bins <= cfg.maxBinSequence;
Raccumul  = Raccumul(isGoodBin,:);
Laccumul  = Laccumul(isGoodBin,:);
RvisGuide = RvisGuide(isGoodBin,:);
LvisGuide = LvisGuide(isGoodBin,:);

seq.bins  = bins(isGoodBin);

%% normalize to max for display
Raccumul  = Raccumul ./repmat(max(abs(Raccumul)),[sum(isGoodBin) 1]);
Laccumul  = Laccumul ./repmat(max(abs(Laccumul)),[sum(isGoodBin) 1]);
RvisGuide = RvisGuide ./repmat(max(abs(RvisGuide)),[sum(isGoodBin) 1]);
LvisGuide = LvisGuide ./repmat(max(abs(LvisGuide)),[sum(isGoodBin) 1]);

%% COM
seq.accumul.R.COM            = COM(Raccumul',seq.bins);
seq.accumul.L.COM            = COM(Laccumul',seq.bins);
seq.visGuide.R.COM           = COM(RvisGuide',seq.bins);
seq.visGuide.L.COM           = COM(LvisGuide',seq.bins);

%% sort COM
[~,seq.accumul.R.sortIdx]    = sort(seq.accumul.R.COM);
[~,seq.accumul.L.sortIdx]    = sort(seq.accumul.L.COM);
[~,seq.visGuide.R.sortIdx]   = sort(seq.visGuide.R.COM);
[~,seq.visGuide.L.sortIdx]   = sort(seq.visGuide.L.COM);

%% sort dff for display

% by side
seq.accumul.R.sortedbyR.dff  = Raccumul(:,seq.accumul.R.sortIdx);
seq.accumul.R.sortedbyL.dff  = Raccumul(:,seq.accumul.L.sortIdx);
seq.accumul.L.sortedbyR.dff  = Laccumul(:,seq.accumul.R.sortIdx);
seq.accumul.L.sortedbyL.dff  = Laccumul(:,seq.accumul.L.sortIdx);

seq.visGuide.R.sortedbyR.dff = RvisGuide(:,seq.visGuide.R.sortIdx);
seq.visGuide.R.sortedbyL.dff = RvisGuide(:,seq.visGuide.L.sortIdx);
seq.visGuide.L.sortedbyR.dff = LvisGuide(:,seq.visGuide.R.sortIdx);
seq.visGuide.L.sortedbyL.dff = LvisGuide(:,seq.visGuide.L.sortIdx);

% by maze
seq.visGuide.R.sortedbyRaccumul.dff = RvisGuide(:,seq.accumul.R.sortIdx);
seq.visGuide.L.sortedbyLaccumul.dff = LvisGuide(:,seq.accumul.L.sortIdx);

%% summarize: correlations between sorted indices and avg COM, as well as diff for stats

% corr between sort indices 
seq.accumul.SORTcorr_RvsL      = corr(seq.accumul.R.sortIdx,seq.accumul.L.sortIdx);
seq.visGuide.SORTcorr_RvsL     = corr(seq.visGuide.R.sortIdx,seq.visGuide.L.sortIdx);
seq.SORTcorr_accumulVSvisGuide = corr([seq.visGuide.R.sortIdx; seq.visGuide.L.sortIdx],[seq.accumul.R.sortIdx; seq.accumul.L.sortIdx]);

% pval for difference across individual sorts
seq.accumul.SORTp_RvsL         = signrank(seq.accumul.R.sortIdx,seq.accumul.L.sortIdx);
seq.visGuide.SORTp_RvsL        = signrank(seq.visGuide.R.sortIdx,seq.visGuide.L.sortIdx);
seq.SORTp_accumulVSvisGuide    = signrank([seq.visGuide.R.sortIdx; seq.visGuide.L.sortIdx],[seq.accumul.R.sortIdx; seq.accumul.L.sortIdx]);

% corr between unsorted COMs 
seq.accumul.COMcorr_RvsL       = corr(seq.accumul.R.COM,seq.accumul.L.COM);
seq.visGuide.COMcorr_RvsL      = corr(seq.visGuide.R.COM,seq.visGuide.L.COM);
seq.COMcorr_accumulVSvisGuide  = corr([seq.visGuide.R.COM; seq.visGuide.L.COM],[seq.accumul.R.COM; seq.accumul.L.COM]);

% corr between overall activity
seq.accumul.OVERALLcorr_RvsL       = corr(Raccumul(:),Laccumul(:));
seq.visGuide.OVERALLcorr_RvsL      = corr(RvisGuide(:),LvisGuide(:));
seq.OVERALLcorr_accumulVSvisGuide  = corr([Raccumul(:); Laccumul(:)],[RvisGuide(:); LvisGuide(:)]);

% pval for difference across individual pxl / rois
seq.accumul.COMp_RvsL          = signrank(seq.accumul.R.COM,seq.accumul.L.COM);
seq.visGuide.COMp_RvsL         = signrank(seq.visGuide.R.COM,seq.visGuide.L.COM);
seq.COMp_accumulVSvisGuide     = signrank([seq.visGuide.R.COM; seq.visGuide.L.COM],[seq.accumul.R.COM; seq.accumul.L.COM]);


%% plot
if plotFlag
  fh   = figure;
  wf   = widefieldParams;
  nROI = size(seq.accumul.R.sortedbyR.dff,2);
  wf.applyFigDefaults(fh,[4 2],'k')
  
  subplot(2,3,1)
  imagesc(seq.bins,1:nROI,seq.accumul.R.sortedbyR.dff',[-1 1])
  wf.applyAxisDefaults(gca,'w'); 
  wf.applyAxisLbls(gca,'y (cm)','ROI','accumul R sort by R')
  
  subplot(2,3,2)
  imagesc(seq.bins,1:nROI,seq.accumul.R.sortedbyL.dff',[-1 1])
  wf.applyAxisDefaults(gca,'w'); 
  wf.applyAxisLbls(gca,'y (cm)','ROI','accumul R sort by L')
  
  subplot(2,3,4)
  imagesc(seq.bins,1:nROI,seq.visGuide.R.sortedbyR.dff',[-1 1])
  wf.applyAxisDefaults(gca,'w'); 
  wf.applyAxisLbls(gca,'y (cm)','ROI','vis guide R sort by R')
  
  subplot(2,3,5)
  imagesc(seq.bins,1:nROI,seq.visGuide.R.sortedbyL.dff',[-1 1])
  wf.applyAxisDefaults(gca,'w'); 
  wf.applyAxisLbls(gca,'y (cm)','ROI','vis guide R sort by L')
  
  subplot(2,3,6)
  imagesc(seq.bins,1:nROI,seq.visGuide.R.sortedbyRaccumul.dff',[-1 1])
  wf.applyAxisDefaults(gca,'w'); 
  wf.applyAxisLbls(gca,'y (cm)','ROI','vis guide R sort by R acc.')
  
  subplot(2,3,3)
  set(gca,'xtick',[],'ytick',[],'color','k')
  xlim([0 1]); ylim([0 1])
  text(.05,.5,sprintf('COM corr (accumul, R vs L) = %1.2f\nCOM corr (visGuide, R vs L) = %1.2f\nCOM corr (accumul vs visGuide) = %1.2f',...
       seq.accumul.COMcorr_RvsL,seq.visGuide.COMcorr_RvsL,seq.COMcorr_accumulVSvisGuide),...
       'verticalAlignment','middle','color','w','fontsize',11)
end