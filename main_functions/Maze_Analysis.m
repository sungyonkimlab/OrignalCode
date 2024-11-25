%% O-Maze <3 or Noldus video data w/ LED sync analysis
clear all;close all

%load data
fpObj = FPObjMake;
if   fpObj(1).loaded == 0
    
    %data pre processing
    guiOut = fpGUI_2;
    fpObj = applyParameters(fpObj,guiOut);
    
    
    fpObj = getTTLOnOffTime(fpObj);
    fpObj = getEventWindowIdx(fpObj);
    trimGuiOut = trimmingGUI_2;
    %trimming data and get dFF
    fpObj = setDataTrimming(fpObj,trimGuiOut);
    fpObj = getTimeVectors(fpObj);
    fpObj = applyTrimmingOffset(fpObj);
    fpObj = calculatedFF_choice(fpObj);
    
    %trimmed data via TTL Onset
    fpObj = trimdFF_TTLOnset(fpObj);
    
    %% Noldus data Analysis
    %Noldus Data Loading
    %loadNoldusFile(fpObj, Meta data column number for noldus file path)
    tic
    fpObj = getVideoOnOffIdx(fpObj);
    toc
    fpObj = loadNoldusFile(fpObj,6);
    saveFPObj(fpObj)

end

%Analyze and plot 2D heatmap with mouse position
fpObj = getNoldusPosition(fpObj);
fpObj = interpolatingNoldusPosition(fpObj);
% fpObj = calculate2dMeandFF(fpObj);
fpObj = calculateMeanandSTDofdFF(fpObj);
fpObj = calculate2dNormdFF(fpObj);

%Plot 2d heatmap
% plot2dHeatmap(fpObj);

% Light-Dark Box Analysis
% fpObj = getZoneTime(fpObj,zone number,vector of colum number about zone entrance);
% Ex) fpObj = fpObj = getZoneTime(fpObj,2,[15 16]);    Column 15, 16 is
% inZone and inZone2
% fpObj = getZoneTime(fpObj,2,[10 11]);
% fpObj = getFPDataIdx(fpObj);
