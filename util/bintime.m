function [binmat,bintaxis] = bintime(mat,taxis,binFactor,meanOrSum)

% [binmat,bintaxis] = bintime(mat,taxis,binFactor,meanOrSum)
% mat is trials x time 
% taxis is time stamps corresponding to cols in mat
% bin factor is factor to bin by (default 4)
% meanOrSum (default 'mean'), to average or sum over timepoints

%%
if nargin < 2; taxis     = [];      end
if nargin < 3; binFactor = 4;       end
if nargin < 4; meanOrSum = 'mean';  end

%%
if binFactor == 1
    binmat   = mat;
    bintaxis = taxis;
else
    bins  = 1:binFactor:size(mat,2);
    if bins(end) < size(mat,2); bins(end+1) = size(mat,2); end
    
    binmat = zeros(size(mat,1),length(bins)-1);
    if ~isempty(taxis); bintaxis = zeros(1,length(bins)-1); end
    
    for bb = 1:length(bins)-1
      switch meanOrSum
        case 'mean'
          binmat(:,bb) = nanmean(mat(:,bins(bb):bins(bb+1)),2);
          
        case 'sum'
          binmat(:,bb) = nansum(mat(:,bins(bb):bins(bb+1)),2);
          
      end
      if ~isempty(taxis)
        bintaxis(bb) = mean(taxis(:,bins(bb):bins(bb+1)));
      end
    end
end
    