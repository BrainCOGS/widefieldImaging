classdef widefield_recLs

  properties (Constant)
    
    % list of gcamp6f recordings during behavior, with isosbestic
    % correction
    BVrecs_gcamp = {          ...
              'ai2/20170125'; ...
              'ai2/20170126'; ...
              'ai3/20170126'; ...
              'ai3/20170201'; ...
              'ai3/20170202'; ...
              'ai3/20170224'; ...
              'ai3/20170227'; ...
              'ai5/20170224'; ...
              'ai5/20170228'; ...
              'ai5/20170313'; ...
              'ai7/20170808'; ...
              'ai7/20170811'; ...
              'ai7/20170814'; ...
              'ai7/20170914'; ...
              'ai7/20171219'; ...
              'ai9/20180327'; ...
              'ai9/20180328'; ...
              'ai9/20180329'; ...
              'ai9/20180404'; ...
              'ai9/20180406'; ...
              'ai10/20180403'; ...
              'ai10/20180404'; ...
              'ai10/20180405'; ...
              'ai10/20180406'; ...
              'ai10/20180409'; ...
              };
    
    % list of gcamp6f recordings during running in the dark, with isosbestic
    % correction        
    BVrecs_spont = {                 ...
              'ai2/20170410_spont';  ...
              'ai3/20170223_spont';  ...
              'ai5/20170223_spont';  ...
              'ai7/20170810_spont';  ...
              'ai9/20180404_spont';  ...
              'ai10/20180404_spont'; ...
              }
    
    % list of yfp recordings during behavior, with isosbestic
    % correction        
    BVrecs_yfp = {...
              'ty1/20170321'; ...
              'ty1/20170331'; ...
              'ty2/20170125'; ...
              'ty2/20170127'; ...
              'ty2/20170131'; ...
              'ty2/20170201'; ...
              'ty2/20170313'; ...
              };
            
    BVrecs   = [widefield_recLs.BVrecs_gcamp; widefield_recLs.BVrecs_yfp];
    BVrecsAll= [widefield_recLs.BVrecs_gcamp; widefield_recLs.BVrecs_yfp; widefield_recLs.BVrecs_spont];
    
    % list of gcamp6f recordings during visual mapping
    visMaps  = { 'ai2/20160721_visMap/retinotopy_20160818_1726/ai2_20160721_visMap_bilateral_000__20160818_1726_retinotopy.mat';             ...
                 'ai3/20160803_visMap/retinotopy_20170630_1228/ai3_ai3_20160803_visMap_bilateral_000__20170630_1228_retinotopy.mat';         ...
                 'ai5/20170411_visMap_bilateral/retinotopy_20170628_1306/ai5_20170411_visMap_bilateral_000__20170628_1306_retinotopy.mat';   ...
                 'ai7/20170814_visMap_bilateral/retinotopy_20170815_1547/ai7_visMapping_bilateral_000__20170815_1547_retinotopy.mat';        ...
                 'ai9/20180404_visMap_bilateral/retinotopy_20180406_1344/ai9_20180404_visMap_bilateral_000__20180406_1344_retinotopy.mat';   ...
                 'ai10/20180404_visMap_bilateral/retinotopy_20180406_1417/ai10_20180404_visMap_bilateral_000__20180406_1417_retinotopy.mat'; ...
                 };
    visMaps_isBilateral = [true; true; true; true; true; true];
    
    % list of gcamp6f recordings during whisker pad mapping
    airpuffMaps = { 'ai2/20170410_whiskersR/map_ai2_whiskersR_20170410_1330.mat';   ...
                    'ai3/20170411_whiskersR/map_ai3_whiskersR_20170411_1340.mat';   ...
                    'ai5/20170411_whiskersR/map_ai5_whiskersR_20170411_1141.mat';   ...
                    'ai7/20170814_whiskersR/map_ai7_whiskersR_20170814_1555.mat';   ...
                    'ai9/20180404_whiskersR/map_ai9_whiskersR_20180404_1454.mat';   ...
                    'ai10/20180404_whiskersR/map_ai10_whiskersR_20180404_1625.mat'; ...
                  }
    airpuffMaps_isBilateral = [true; true; true; true; true; true]
    airpuffMaps_stimWhisker = ['R';'R';'R';'R';'R';'R']      
  end
end