function [recls,rootdir] = getFullRecPath(recls)

% [recls,rootdir] = getFullRecPath(recls)
% recls can be a cell array with partial or full paths, or just the mouse
% name, in which case it will load from widefield_recLs

if ischar(recls)
  fullrecls = widefield_recLs.BVrecs;
  isTarget  = cellfun(@(x)(contains(x,recls)),fullrecls);
  recls     = fullrecls(isTarget);
end

if ~isFullPath(recls{1})
  wf        = widefieldParams;
  rootdir   = wf.getRootDir(isThisSpock);
  recls     = cellfun(@(x)(formatFilePath([rootdir x])),recls,'Uniformoutput',false);
else
  mouseID   = mouseAndDateFromFileName(recls{1});
  rootdir   = recls{1}(1:strfind(recls{1},mouseID)-1);
end