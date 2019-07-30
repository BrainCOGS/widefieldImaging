function strout = removeUnderscores(strin)

% strout = removeUnderscores(strin)
% replaces underscores by spaces in string

idx         = strfind(strin,'_');
strout      = strin;
strout(idx) = ' ';