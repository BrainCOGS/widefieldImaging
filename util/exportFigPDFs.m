function exportFigPDFs(fn)

% exportFigPDFs(fn)
% exports all open figures to a pdf (filename fn, with .pdf extension) and
% closes them

if ~isempty(dir(fn)); delete(fn); end

figls = get(0,'children'); % print to pdf
for ii = length(figls):-1:1
  figure(figls(ii))
  export_fig(fn,'-q101','-append')
end
close all