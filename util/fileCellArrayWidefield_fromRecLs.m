function flp = fileCellArrayWidefield_fromRecLs(recls,cohort,spockFlag)

% [flp,fl,mouseID] = fileCellArrayWidefield(mice,cohort)
% generates cell array of file paths for desired recs
% INPUT
%   recls: cell array of rec names 
%   cohort: name of mouse cohort, def: blocksReboot
%
% OUTPUT
%   flp and fl are structures containing full paths/just file names for
%       .behav: virmen log files
%       .widefield:   lsr log files
%   mouseID mouse to which each entry in the structure fields belongs
%
% LP june 2017

%% defaults
if nargin < 2 || isempty(cohort)
  cohort    = 'blocksReboot';
end
if nargin < 3
  spockFlag = true;
end

fprintf('generating file list...\n')

%% get paths
if spockFlag
  rootdir_wf    = widefieldParams.serverpath_spock;
  rootdir_bhv   = [widefieldParams.serverpathbehav_spock cohort '\data\'];
else
  if ispc
    rootdir_wf  = widefieldParams.serverpath_pc;
    rootdir_bhv = [widefieldParams.serverpathbehav_pc cohort '\data\'];
  else
    rootdir_wf  = widefieldParams.serverpath_mac;
    rootdir_bhv = [widefieldParams.serverpathbehav_mac cohort '\data\'];
  end
end

%% loop thorugh recs
flp = cell(numel(recls),1);

for iRec = 1:numel(recls)
  
  % info file path (for syncing)
  flp{iRec}.widefield = formatFilePath([rootdir_wf recls{iRec} '/info.mat'],false);
  [rmouse,rdate]      = mouseAndDateFromFileName((recls{iRec}));
  
  % virmen log paths
  cd(formatFilePath([rootdir_bhv rmouse]))
  subfl           = dir('*WideFieldCohort*mat');
  subfl           = {subfl(:).name}';
  recidx          = cellfun(@(x)(~isempty(strfind(x,rdate))),subfl);
  flp{iRec}.behav = formatFilePath([rootdir_bhv rmouse '/' subfl{recidx}],false);
  
end

