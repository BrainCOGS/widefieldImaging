rd = '/jukebox/braininit/RigData/VRwidefield/widefield/';
recls = widefield_recLs.BVrecs_gcamp;
model_names = {'original','+corr','+auto_regr','+both'};

cfg_all{1}.predList      = {'tow_L','tow_R','\Delta','y','\theta','d\theta/dt','speed','ch','prevch','prevrw'};
cfg_all{1}.predLagCm     = {[0 100],[0 100],[0 50],[0 0],[0 0],[30 30],[30 30],[0 0],[0 0],[0 0]};                       

cfg_all{2}.predList      = {'tow_L','tow_R','\Delta','y','\theta','d\theta/dt','speed','ch','prevch','prevrw','ROI'};
cfg_all{2}.predLagCm     = {[0 100],[0 100],[0 50],[0 0],[0 0],[30 30],[30 30],[0 0],[0 0],[0 0],[0 0]};                       

cfg_all{3}.predList      = {'tow_L','tow_R','\Delta','y','\theta','d\theta/dt','speed','ch','prevch','prevrw','auto'};
cfg_all{3}.predLagCm     = {[0 100],[0 100],[0 50],[0 0],[0 0],[30 30],[30 30],[0 0],[0 0],[0 0],[0 100]};                       

cfg_all{4}.predList      = {'tow_L','tow_R','\Delta','y','\theta','d\theta/dt','speed','ch','prevch','prevrw','ROI','auto'};
cfg_all{4}.predLagCm     = {[0 100],[0 100],[0 50],[0 0],[0 0],[30 30],[30 30],[0 0],[0 0],[0 0],[0 0],[0 100]};                       


nROI     = 16;
nRecs    = numel(recls);
nModels  = numel(model_names);
accuracy = zeros(nRecs,nROI,nModels);

for iRec = 1:nRecs
  cd([rd recls{iRec}])
  for iModel = 1:nModels
    fprintf('\n\nREC %02d/%02d, MODEL %d/%d\n\n',iRec,nRecs,iModel,nModels)
    if iModel > 1; saveFlag = true; else; saveFlag = false; end
    this_glm = dffGLM(pwd,[],1,0,cfg_all{iModel},{},saveFlag,0);
    accuracy(iRec,:,iModel) = this_glm.accuracy;
  end
end

clear iRec iModel this_glm
cd(rd)
save dffGLM_autoregr_model_comp