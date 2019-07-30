function [fl,flp] = generateFilelist(rec)

% [fl,flp] = generateFilelist(rec)
% generates cell array with all .tif file names (fl) or full paths(flp)

% make sure it's full path
rec = formatFilePath(rec);
if isdir(rec)
  fp = rec;
else
  if ispc
    fp = sprintf('%s%s',widefieldParams.savepath_pc,rec);
  else
    fp = sprintf('%s%s',widefieldParams.savepath_mac,rec);
  end
end

% list
fl  = dir([fp '*tif']);
fl  = {fl.name};
flp = cellfun(@(x)([fp x]),fl,'UniformOutput', false);
