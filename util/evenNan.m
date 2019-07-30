% imstack = evenNan(imstack)
% make nans uniform across image stack (useful for aligining file chunks 
% with different vasculature / motion correction). If a pixel contains nans
% it will be nan'd across the stack
function imstack = evenNan(imstack)

[nX,nY,nZ] = size(imstack);
hasNan     = false(nX,nY);

for iRow = 1:nX
  for iCol = 1:nY
    hasNan(iRow,iCol) = any(isnan(imstack(iRow,iCol,:)));
  end
end

imstack(repmat(hasNan,[1 1 nZ])) = nan;

end