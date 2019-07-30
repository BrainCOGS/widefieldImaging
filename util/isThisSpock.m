function isSpock = isThisSpock

% isSpock = isThisSpock
% returns true if current envirnment is within PNI computing cluser

if (contains(pwd,'smb')            || ...
    contains(pwd,'usr/people')     || ...
    contains(pwd,'mnt')            || ...
    contains(pwd,'jukebox'))          ...
    && ~ispc
    
  isSpock = true;
  
else
  
  isSpock = false;
  
end
