clear all; close all;

%load data
fpObj = FPObjMake;
%data pre processing
if  fpObj(1).loaded == 0
    
    guiOut = fpGUI_2;
    fpObj = applyParameters(fpObj,guiOut);
    fpObj = getTTLOnOffTime(fpObj);
    trimGuiOut = trimmingGUI_2;

    %trimming data and get dFF
    fpObj = setDataTrimming(fpObj,trimGuiOut);
    fpObj = getEventWindowIdx(fpObj);
    fpObj = getTimeVectors(fpObj);
    fpObj = applyTrimmingOffset(fpObj);
    fpObj = calculatedFF_choice(fpObj);
    
    %visualize
    fpObj = calculateFTT_dFF(fpObj);
    saveFPObj(fpObj)
end


plotFTTHeatmap(fpObj);

