function newcoord = reflectCoordOverMidline(coord,midline,toSide,imSize)

% newcoord = reflectCoordOverMidline(coord,midline,toSide,imSize)
% find the homotopic coordinate given the current midline coordinates

if nargin < 4
  imSize = [128 128];
end

a        = midline.a;
c        = midline.b;
b        = -1;

newcoord = nan(size(coord));
for iPxl = 1:size(coord,1)
  if isfinite(a) && isfinite(c)
    x    = coord(iPxl,2);
    y    = coord(iPxl,1); 

    % distance between this pixel and the midline
    dist = abs(a*x + b*y + c) / sqrt(a^2 + b^2);

    % point on midline closest to (x,y)
    x0   = (b*(b*x - a*y) -a*c) / (a^2 + b^2);
    y0   = (a*(-b*x + a*y)-b*c) / (a^2 + b^2);

    % slope of line that connects (x,y) and (x0,y0)
    if strcmpi(toSide,'R')
      m  = (y0-y)/(x0-x);
    else
      m  = (y-y0)/(x-x0);
    end
    int  = y0 - m*x0;

    % now find point on this line that lies at distance dist from midline
    if strcmpi(toSide,'R')
      xn = x0 + dist / sqrt(1+m^2);
      yn = m*xn + int;
    else
      xn = x0 - dist / sqrt(1+m^2);
      yn = int + m*xn;
    end
    
  else % if a /c are inf it means midline is parallel to y
    yn   = coord(iPxl,1); 
    xmid = midline.p(1,1);
    dist = xmid - coord(iPxl,2);
    xn   = xmid + dist;
  end
  
  newcoord(iPxl,:) = [yn xn];
end

newcoord = round(newcoord); 
delidx   = newcoord(:,1) < 1 | newcoord(:,1) > imSize(1) | ...
           newcoord(:,2) < 1 | newcoord(:,2) > imSize(2) | ...
           isnan(newcoord(:,1)) | isnan(newcoord(:,2));
newcoord = newcoord(~delidx,:);