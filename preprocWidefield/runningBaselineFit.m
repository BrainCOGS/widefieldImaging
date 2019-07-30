function [f0,cf] = runningBaselineFit(B,V,ws,frameRate,method,doSmooth)

% [f0,cf] = runningBaselineFit(B,V,ws,frameRate,method,doSmooth)
% calculates pixel-wise heuristic correct factor cf as the slope of the
% linear regression B_hat = F0 = cf*V + a and returns a corrected F0 for the blue
% channel. 
% B is blue fluorescence , V is violet fluorescence, frame rate is
% recording frame rate, recommended method is 'overall', doSmooth true to
% first smooth V
% called by dffFrom2Channels()

%% defaults
if nargin < 3
  ws = 30; %sec
end
if nargin < 4
  try
    load info frameRate
  catch
    frameRate = 20;
  end
end
if nargin < 5
  method = 'overall'; 
end
if nargin < 6
  doSmooth = 1;
end

%% fit baseline
tic
fprintf('fitting baseline... ')

fr      = frameRate / 2;
f0      = nan(size(B));
nframes = size(V,1);
npxls   = size(V,2);
cf      = cell(1,npxls);

if doSmooth
  vi  = V;
  gws = round(.4*fr);
  V   = approxGaussianBlur1D(vi,gws,1);
  V([1:gws*3 nframes-gws*3:nframes],:) = V([1:gws*3 nframes-gws*3:nframes],:);
end

switch method
  case 'overall'
    nframesB    = size(B,1);
    nfit        = min([nframes nframesB]);
    for iPxl = 1:npxls
      cf{iPxl}    = fit(V(100:nfit-100,iPxl),B(100:nfit-100,iPxl),'poly1');
      f0(:,iPxl)  = fiteval(cf{iPxl},V(:,iPxl));
    end
  
  case 'offset'
    nframesB    = size(B,1);
    nfit        = min([nframes nframesB]);
    for iPxl = 1:npxls
      cf{iPxl}    = fit(V(100:nfit-100,iPxl),B(100:nfit-100,iPxl),fittype('x+p2'),'startPoint',-7e4);
      f0(:,iPxl)  = fiteval(cf{iPxl},V(:,iPxl));
    end
    
  case 'chunks'
    bins = 1:round(ws*fr):nframes;
    if bins(end) < nframes; bins(end) = nframes; end
    for iW = 1:numel(bins)-1
      v                        = V(bins(iW):bins(iW+1),:);
      b                        = B(bins(iW):bins(iW+1),:);
      for iPxl = 1:npxls
        cf{iPxl}{iW}                  = fit(v(:,iPxl),b(:,iPxl),'poly1');
        f0(bins(iW):bins(iW+1),iPxl)  = fiteval(cf{iPxl},v(:,iPxl));
      end
    end
  case 'sliding'
    for iW = ceil(fr*ws/2):nframes-round(fr*ws/2)
      idx              = iW-floor(fr*ws/2)+1:iW+floor(fr*ws/2)-1;
      v                = V(idx,:);
      b                = B(idx,:);
      
      for iPxl = 1:npxls
        cf{iPxl}{iW}     = fit(v(:,iPxl),b(:,iPxl),'poly1');
        f0(iW,iPxl)      = fiteval(cf{iPxl},V(iW,iPxl));
        
        if iW == ceil(fr*ws/2)
          v                     = V(1:iW-1,iPxl);
          f0(1:iW-1,iPxl)       = fiteval(cf{iPxl} ,v);
        elseif iW == numel(V)-ceil(fr*ws/2)
          v                     = V(iW+1:end,iPxl);
          f0(iW+1:end,iPxl)     = fiteval(cf{iPxl} ,v);
        end
      end
    end
end

fprintf('done after %1.1f min\n',toc/60)
