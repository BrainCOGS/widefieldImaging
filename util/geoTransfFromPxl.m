function tform = geoTransfFromPxl(shiftVec)

% tform = geoTransfFromPxl(shiftVec)
% create affin2d transfomation objet from a vector (shiftVec) of known image
% translation, rotation and scale, 
% [colshift rowshift angle scale], typically from xRecAlign_manual
% LP July 2017

xshift = shiftVec(2);
yshift = shiftVec(1);
scale  = shiftVec(4);
ang    = deg2rad(shiftVec(3));

Tmat   = [scale      0        0;     ...
          0         scale     0;     ...
          xshift    yshift    1];
       
Rmat   = [cos(ang)  sin(ang)  0;     ...
          -sin(ang) cos(ang)  0;     ...
          0         0         1];
        
tform  = affine2d(Tmat*Rmat);