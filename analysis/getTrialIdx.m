function trialidx = getTrialIdx(logSumm,trialtype,maze,perfTh,combineVisGuide,includeAll)

% trialidx = getTrialIdx(logSumm,trialtype,maze,perfTh,combineVisGuide,includeAll)
% returns indices for trials meeting criteria in inputs
%   logSumm is falttened behavior log
%   trialtype is string, any combination of
%        [correct-error][easy-hard][R-L][_distractor][_nodistractor]. Refer
%        to body of function for more details
%   maze is maze number, maze < 0 will select all mazes
%   perfTh: blockwise performance threshold [0 1] to include trial blocks (default 0.6)
%   combineVisGuide to combine visually-guided trials with and without distractors (default true)
%   includeAll: true to relax all bad trial exclusion criteria (default false)

%%
if nargin < 3 || isempty(maze)
  maze = [];
end
if nargin < 4 || isempty(perfTh)
  perfTh = widefieldParams.perfTh;
end
if nargin < 5 || isempty(combineVisGuide)
  combineVisGuide = true; % visual guide with and without distractors?
end
if nargin < 6 || isempty(includeAll)
  includeAll      = false; % relax bad trial criteria?
end

if combineVisGuide %&& ~isempty(strfind(char(logSumm.info.protocol),'Condensed'))
  logSumm.currMaze(logSumm.currMaze == 12) = 4;
end
if isempty(maze)
  maze = max(logSumm.currMaze);
end

%%
trialDifficulty = 1 - (abs(logSumm.nCues_RminusL) ./ logSumm.nCues_total);


switch trialtype
  case 'all'
    idx = 1:numel(logSumm.choice);
  case 'correct'
    idx = find(logSumm.trialType == logSumm.choice);
  case 'error'
    idx = find(logSumm.trialType ~= logSumm.choice);
  case 'R'
    idx = find(logSumm.trialType == analysisParams.rightCode);
  case 'L'
    idx = find(logSumm.trialType == analysisParams.leftCode);
  case 'correctR'
    idx = find(logSumm.trialType == logSumm.choice & logSumm.choice == analysisParams.rightCode);
  case 'correctL'
    idx = find(logSumm.trialType == logSumm.choice & logSumm.choice == analysisParams.leftCode);
  case 'errorR'
    idx = find(logSumm.trialType ~= logSumm.choice & logSumm.choice == analysisParams.leftCode);
  case 'errorL'
    idx = find(logSumm.trialType ~= logSumm.choice & logSumm.choice == analysisParams.rightCode);
  case 'correct_nodistractor'
    idx = find(logSumm.trialType == logSumm.choice & ((logSumm.trialType == analysisParams.rightCode & logSumm.nCues_L == 0) | (logSumm.trialType == analysisParams.leftCode & logSumm.nCues_R == 0)));
  case 'error_nodistractor'
    idx = find(logSumm.trialType ~= logSumm.choice & ((logSumm.trialType == analysisParams.rightCode & logSumm.nCues_L == 0) | (logSumm.trialType == analysisParams.leftCode & logSumm.nCues_R == 0)));
  case 'correct_distractor'
    idx = find(logSumm.trialType == logSumm.choice & logSumm.nCues_R > 0 & logSumm.nCues_L > 0);
  case 'error_distractor'
    idx = find(logSumm.trialType ~= logSumm.choice & logSumm.nCues_R > 0 & logSumm.nCues_L > 0);
  case 'correctR_distractor'
    idx = find(logSumm.trialType == logSumm.choice & logSumm.choice == analysisParams.rightCode & logSumm.nCues_L > 0);
  case 'correctL_distractor'
    idx = find(logSumm.trialType == logSumm.choice & logSumm.choice == analysisParams.leftCode & logSumm.nCues_R > 0);
  case 'errorR_distractor'
    idx = find(logSumm.trialType ~= logSumm.choice & logSumm.choice == analysisParams.leftCode & logSumm.nCues_L > 0);
  case 'errorL_distractor'
    idx = find(logSumm.trialType ~= logSumm.choice & logSumm.choice == analysisParams.rightCode & logSumm.nCues_R > 0);
  case 'correctR_nodistractor'
    idx = find(logSumm.trialType == logSumm.choice & logSumm.choice == analysisParams.rightCode & logSumm.nCues_L == 0);
  case 'correctL_nodistractor'
    idx = find(logSumm.trialType == logSumm.choice & logSumm.choice == analysisParams.leftCode & logSumm.nCues_R == 0);
  case 'errorR_nodistractor'
    idx = find(logSumm.trialType ~= logSumm.choice & logSumm.choice == analysisParams.leftCode & logSumm.nCues_L == 0);
  case 'errorL_nodistractor'
    idx = find(logSumm.trialType ~= logSumm.choice & logSumm.choice == analysisParams.rightCode & logSumm.nCues_R == 0);
  case 'easy'
%     idx = find(trialDifficulty < nanmedian(trialDifficulty)); % old definition 
    idx = find(trialDifficulty < nanmedian(trialDifficulty));
  case 'hard'
%     idx = find(trialDifficulty > nanmedian(trialDifficulty));
    idx = find(trialDifficulty > nanmedian(trialDifficulty));
  case 'easyL'
     idx = find(logSumm.trialType == analysisParams.leftCode & trialDifficulty < nanmedian(trialDifficulty));
  case 'hardL'
     idx = find(logSumm.trialType == analysisParams.leftCode & trialDifficulty > nanmedian(trialDifficulty));
  case 'easyR'
     idx = find(logSumm.trialType == analysisParams.rightCode & trialDifficulty < nanmedian(trialDifficulty));
  case 'hardR'
     idx = find(logSumm.trialType == analysisParams.leftCode & trialDifficulty > nanmedian(trialDifficulty));
  case 'correct_easy'
    idx = find(logSumm.trialType == logSumm.choice & trialDifficulty < nanmedian(trialDifficulty));
  case 'error_easy'
    idx = find(logSumm.trialType ~= logSumm.choice & trialDifficulty < nanmedian(trialDifficulty));
  case 'correct_hard'
    idx = find(logSumm.trialType == logSumm.choice & trialDifficulty > nanmedian(trialDifficulty));
  case 'error_hard'
    idx = find(logSumm.trialType ~= logSumm.choice & trialDifficulty > nanmedian(trialDifficulty));
  case 'correctL_easy'
    idx = find(logSumm.trialType == logSumm.choice & logSumm.trialType == analysisParams.leftCode & trialDifficulty < nanmedian(trialDifficulty));
  case 'errorL_easy'
    idx = find(logSumm.trialType ~= logSumm.choice & logSumm.trialType == analysisParams.leftCode & trialDifficulty < nanmedian(trialDifficulty));
  case 'correctL_hard'
    idx = find(logSumm.trialType == logSumm.choice & logSumm.trialType == analysisParams.leftCode & trialDifficulty > nanmedian(trialDifficulty));
  case 'errorL_hard'
    idx = find(logSumm.trialType ~= logSumm.choice & logSumm.trialType == analysisParams.leftCode & trialDifficulty > nanmedian(trialDifficulty));
  case 'correctR_easy'
    idx = find(logSumm.trialType == logSumm.choice & logSumm.trialType == analysisParams.rightCode & trialDifficulty < nanmedian(trialDifficulty));
  case 'errorR_easy'
    idx = find(logSumm.trialType ~= logSumm.choice & logSumm.trialType == analysisParams.rightCode & trialDifficulty < nanmedian(trialDifficulty));
  case 'correctR_hard'
    idx = find(logSumm.trialType == logSumm.choice & logSumm.trialType == analysisParams.rightCode & trialDifficulty > nanmedian(trialDifficulty));
  case 'errorR_hard'
    idx = find(logSumm.trialType ~= logSumm.choice & logSumm.trialType == analysisParams.rightCode & trialDifficulty > nanmedian(trialDifficulty));
end

if includeAll
  goodtrials   = 1:numel(logSumm.currMaze);
else
  goodtrials   = find(~isBadTrial(logSumm,[],false,perfTh,true));
end

if numel(maze) == 1
  if maze < 0
    mazeidx    = 1:numel(logSumm.currMaze);
  else
    mazeidx    = find(logSumm.currMaze == maze);
  end
else
  mazeidx = [];
  for iMz = 1:numel(maze)
    mazeidx = [mazeidx find(logSumm.currMaze == maze(iMz))];
  end
end

trialidx     = intersect(mazeidx,intersect(idx,goodtrials));