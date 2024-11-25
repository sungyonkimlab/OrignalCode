clear all;close all

%load data if directory/file does not exist
% checkFPObj
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
    saveFPObj(fpObj)
end
%% Analysis part for IP injection data

examRange = [-300 2300]; %sec;
fpObj = preprocess_IP_injectionData(fpObj,examRange); %this function retrim data according to TTL pulses that signal 15min window.
