function [idxr,idxl] = assignPxlToHem(coord,midline)

% [idxr,idxl] = assignPxlToHem(coord,midline)
% given midline object, decide if pxls in coord matrix ([pxls x [rows cols]])
% belong to left or right hemisphere
% idxr and idxl are indices of rows in coord matrix

r   = coord(:,1);
c   = coord(:,2);
urs = unique(r);

% left hemisphere
idxl   = []; %left hem pixels
for jj = 1:length(urs)
  % find midline at this AP coordinate
    ym   = urs(jj);
  if isfinite(midline.a) && isfinite(midline.b)
    xm   = (ym-midline.b)/midline.a;
  else
    xm   = midline.p(1,1); % perfectly straight midline
  end
  idxl = [idxl; find(r == ym & c < xm)];
end

idxr = setdiff(1:numel(c),idxl); % right hemisphere pixels
