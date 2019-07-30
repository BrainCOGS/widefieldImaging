function dffROI = avgDffFromROI(imstack,ROI,fn)

% dffROI = avgDffFromROI(imstack,ROI,fn)
% output is a nROI x frame matrix with dff traces
% inputs:   imstack: 3d matrix of dff
%           ROI: cell array with pxl coordinates for each ROI (from autoROI)
%           fn: saved file name, default dffROI.mat
%
% LP Aug 2016 modified from May 2013 (Neuron paper Dan lab version)

if nargin < 1 || isempty(imstack)
  load dffcorrect dffcorrect
  imstack = dffcorrect;
  clear dffcorrect
end
if nargin < 2 || isempty(ROI)
  load ROIfromRef ROI
end
if nargin < 3 || isempty(fn)
  fn = 'dffROI';
end

tic
fprintf('extracting DFF for ROIs...')

% raw fluorescence over ROI and dff
nZ      = size(imstack,3);
nROI    = length(ROI);
dffROI  = zeros(nROI,nZ);

for i = 1:nROI
  ts = zeros(1,nZ); ct = 1;
  for r = 1:length(ROI{i})
    if ~isnan(imstack(ROI{i}(r,1),ROI{i}(r,2),1))
      ts = ts + double(squeeze(imstack(ROI{i}(r,1),ROI{i}(r,2),:))');
      ct = ct + 1;
    end
  end
  dffROI(i,:) = ts./ct;
end

save(fn,'dffROI')

fprintf(' done after %1.1f minutes\n',toc/60)