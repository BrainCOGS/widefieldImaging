function fp = formatFilePath(fp,addSlash)

% fp = formatFilePath(fp,addSlash)
% formats path fp according to current OS environment, addSlash is flag to
% add slash to end of path (default true)

if nargin == 1
  addSlash = 1;
end

bl = strfind(fp,'\');
fl = strfind(fp,'/');

if ispc
  if ~isempty(fl)
    for ii = 1:length(fl)
      fp(fl(ii)) = '\';
    end
  end
  if addSlash && ~strcmpi(fp(end),'\')
    fp(end+1) = '\';
  end
else
  if ~isempty(bl)
    for ii = 1:length(bl)
      fp(bl(ii)) = '/';
    end
  end
  if addSlash && ~strcmpi(fp(end),'/')
    fp(end+1) = '/';
  end
end