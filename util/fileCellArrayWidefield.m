function [flp,fl,mouseID] = fileCellArrayWidefield(mice,cohort)

% [flp,fl,mouseID] = fileCellArrayWidefield(mice,cohort)
% generates cell array of file paths for desired mice and rigs
% INPUT
%   mice: cell array of mouse names (def: {'wt5';'wt6';'vg3';'vg4'})
%   cohort: name of mouse cohort, def: blocksReboot
%
% OUTPUT
%   flp and fl are structures containing full paths/just file names for
%       .behav: virmen log files
%       .widefield:   lsr log files
%   mouseID mouse to which each entry in the structure fields belongs
%
% LP feb 2016

fprintf('generating file list...\n')

% defaults / load variables
defaults   = {analysisParams.mice                              ;   ...
  'blocksReboot'};

inputnames = {'mice';'cohort';'protocol'};
if nargin < length(defaults)
  for i = nargin+1:length(defaults)
    temp = defaults{i};
    eval([inputnames{i} '=temp;']);
  end
end

% handle input types
if ~iscell(mice)     ;  temp = mice; clear mice; mice{1} = temp;   end
if isempty(cohort)   ;  cohort = 'blocksReboot';                   end
if iscell(cohort)    ;  cohort = cohort{1};                        end

fl.behav     = []; flp.behav     = {};
fl.widefield = []; flp.widefield = {};
mouseID      = [];

for a = 1:length(mice)
  
  if ispc
    thisbehav = sprintf('%s%s\\data\\%s',widefieldParams.serverpathbehav_pc,cohort,mice{a});
    thiswide  = sprintf('%s%s\\',widefieldParams.serverpath_pc,mice{a});
  else
    thisbehav = sprintf('%s%s/data/%s',widefieldParams.serverpathbehav_mac,cohort,mice{a});
    thiswide  = sprintf('%s%s/',widefieldParams.serverpath_mac,mice{a});
  end
  
  % virmen log paths
  cd(thisbehav)
  
  subfl = glob('*WideFieldCohort*mat');
  fl.behav = [fl.behav; subfl];
  for f = 1:length(subfl)
    flp.behav{end+1}=formatFile([thisbehav '/' subfl{f}]);
  end
  
  % mouse
  mouseID(end+1:end+length(subfl),:) = repmat(mice{a},[length(subfl),1]);
  
  % laser log paths
  cd(thiswide)
  for f = 1:size(subfl,1)
    
    thispath = subfl(f,:); thispath = thispath{1};
    [dtS,dtE]= regexp(thispath,'[0-9]{8,}'); % look for date in file names (8 digits)
    thisdate = thispath(dtS:dtE);
    
    % look for associated dff file
    if ~isempty(dir(thisdate))
      cd(thisdate)
      if ~isempty(dir('info.mat'))
        flp.widefield{end+1} = formatFilePath([thiswide thisdate '/info.mat'],false);
      else
        flp.widefield{end+1} = [];
      end
    else
      flp.widefield{end+1}   = [];
    end
    
    cd(thiswide)
  end
  
end

% format
flp.behav     = flp.behav';
flp.widefield = flp.widefield';
mouseID       = char(mouseID);
