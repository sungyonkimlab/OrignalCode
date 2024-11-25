function plotBar(fpObj,inspectRange,saveFigures)
%initialization
totalMouseNum = fpObj.totalMouseNum;
%check directory and make dir if does not exists
% if exist('BarGraph') ~= 7
%     [status, msg, msgID] = mkdir('BarGraph');
%     if status == 1
%         cd('BarGraph');
%     end
% else
%     cd('BarGraph');
% end
meanDffStack = [];
meanNonDffStack = [];
examRange = fpObj.examRange;
examRangeIdx = fpObj.examRangeIdx;
samplingRate = round(fpObj.samplingRate);

for numMouse = 1:totalMouseNum
    %% bar graph
    f1 = figure('Units','inch','Position',[1 1 8 8],'Visible','on');
    %initialize for each mouse
    boutIdx = fpObj.idvData(numMouse).boutIdx;
    totalNumBout = fpObj.idvData(numMouse).totalNumBout;
    dFF = fpObj.idvData(numMouse).dFF;
    meanBoutDff = [];
    hozconcatDff = [];
    hozconcat_nonDff = [];
    for boutNum = 1:totalNumBout
        meanBoutDff(boutNum,1) = mean(dFF(boutIdx(boutNum,1):boutIdx(boutNum,2)));
        steBoutDff(boutNum,1) = std(dFF(boutIdx(boutNum,1):boutIdx(boutNum,2)),0,1)/length(dFF(boutIdx(boutNum,1):boutIdx(boutNum,2)));
        hozconcatDff = [hozconcatDff dFF(boutIdx(boutNum,1):boutIdx(boutNum,2))'];
        hozconcat_nonDff = (sum(dFF) - sum(hozconcatDff)) / (size(dFF,1) - size(hozconcatDff,2));
    end
    
    meanDff = mean(hozconcatDff);
    meanNonDff = mean(hozconcat_nonDff);
    %store meanDff and meanNonDff
    
    meanDffStack = [meanDffStack; meanDff];
    meanNonDffStack = [meanNonDffStack;meanNonDff];
    
    subplot(2,2,1)
    boutBar = bar(1,meanDff,'EdgeColor',[0 0 0],'FaceColor',[0 0 1]);
    hold on
    nonBoutBar = bar(2,meanNonDff,'EdgeColor',[0 0 0],'FaceColor',[0 0 0]);
    
    %Xlabel
    set(gca,'xtick',[])
    Labels = {'Bout', 'NonBout'};
    set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
    %Ylabel
    ylabel('\DeltaF/F (%)');
    %ylim,other stuff
    set(gcf,'Color',[1 1 1])
    set(gca,'linewidth',1.6,'FontSize',13,'FontName','Arial')
    set(gca,'box','off')
    set(gca,'TickDir','out'); % The only other option is 'in'
    
    
    title([fpObj.idvData(numMouse).Description '  bar'],'FontSize',6)
    
%     if saveFigures =='y' || saveFigures =='Y'
%         saveas(gcf,[fpObj.idvData(numMouse).Description '  bar.jpg'])
%         saveas(gcf,[fpObj.idvData(numMouse).Description '  bar.svg'])
%     end
    %     export_fig([fpObj.idvData(numMouse).Description ' bar'],'-eps','-jpg')
    
    
    %%
    
    subplot(2,2,2)
    hold on
    bar(meanBoutDff)
    xlabel('Bout Number');
    ylabel('\DeltaF/F (%)');
    xlim([0 totalNumBout+1])
    %ylim,other stuff
    set(gcf,'Color',[1 1 1])
    set(gca,'linewidth',1.6,'FontSize',13,'FontName','Arial','box','off')
    set(gca,'TickDir','out'); % The only other option is 'in'
    
    title({[fpObj.idvData(numMouse).Description ' bout bar'];...
        ['numBout = ' num2str(totalNumBout)]}, 'FontSize',6)
    hold off
    %     export_fig([fpObj.idvData(numMouse).Description ' bout bar'],'-eps','-jpg')
%     if saveFigures =='y' || saveFigures =='Y'
%         saveas(gcf,[fpObj.idvData(numMouse).Description ' bout bar.jpg'])
%         saveas(gcf,[fpObj.idvData(numMouse).Description ' bout bar.svg'])
%     end
    %% Linear Regression
    x = 1:totalNumBout;
    y = meanBoutDff';
    % get regression coefficient
    p = polyfit(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    r = round(corrcoef(x,y),2);%pearson cor
    
    %standard error
    
    
    %plot
    subplot(2,2,3)
    plot(x,y,'o')
    hold on
    
    plot(x,yfit);
    %graph spec
    xlim([.5 totalNumBout+.5])
    xlabel('Bout Number');
    ylabel('\DeltaF/F (%)');
    
    set(gcf,'Color',[1 1 1])
    set(gca,'linewidth',1.6,'FontSize',13,'FontName','Arial','box','off')
    set(gca,'TickDir','out');
%     dim = [.7 .5 .3 .4];
    str = ['pearson cor = ' num2str(r(1,2))];
%     annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
    title({[fpObj.idvData(numMouse).Description ' correlation'];...
        ['numBout = ' num2str(totalNumBout)];...
        [str]}, 'FontSize',6)
    %     export_fig([fpObj.idvData(numMouse).Description ' bout correlation'],'-eps','-jpg')
%     if saveFigures =='y' || saveFigures =='Y'
%         saveas(gcf,[fpObj.idvData(numMouse).Description ' bout correlation.jpg'])
%         saveas(gcf,[fpObj.idvData(numMouse).Description ' bout correlation.svg'])
%     end
    %% inspect range bar graph
    
    %initialize

    inspectRangeIdx = inspectRange*samplingRate; 
    boutIdx = fpObj.idvData(numMouse).boutIdx;
    totalNumBout = fpObj.idvData(numMouse).totalNumBout;
    inspectRangeBoutIdx = [boutIdx(:,1) boutIdx(:,1)] + repmat(inspectRangeIdx,[size(boutIdx,1) 1]);
    firstBoutDffArray = [];
    for boutNum = 1:totalNumBout
        firstBoutDffArray(boutNum,:) = dFF(inspectRangeBoutIdx(boutNum,1):inspectRangeBoutIdx(boutNum,2),1);
    end
    
    firstHalfIdx = abs(inspectRangeIdx(1));
    beforeLickDff = firstBoutDffArray(:, 1 : firstHalfIdx);
    afterLickDff = firstBoutDffArray(:,firstHalfIdx+1:end);
    
    meanBeforeLickDff = mean(beforeLickDff,2);
    meanAfterLickDff = mean(afterLickDff,2);
    %mean of mean
    mOmBeforeLickDff(numMouse) = mean(meanBeforeLickDff);
    mOmAfterLickDff(numMouse) = mean(meanAfterLickDff);
    %ste of mean
    ste_mOmBeforeLickDff = std(meanBeforeLickDff)/sqrt(size(meanBeforeLickDff,1));
    ste_mOmAfterLickDff = std(meanAfterLickDff)/sqrt(size(meanBeforeLickDff,1));
    
    
    
    
    subplot(2,2,4)

    hold on
    
    b_1 = bar(1,mOmBeforeLickDff(numMouse),'FaceColor',[1 1 1],'LineWidth',1.6);
    b_2 = bar(2,mOmAfterLickDff(numMouse),'FaceColor',[0.5 0.5 0.5],'LineWidth',1.6);
    errorbar(1:2,[mOmBeforeLickDff(numMouse) mOmAfterLickDff(numMouse)],[ste_mOmBeforeLickDff ste_mOmAfterLickDff],'.','Color','k','linewidth',1.3)
    set(gca,'xtick',[])
    Labels = {'Before Lick', 'After Lick'};
    set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
    % set(h(1),'FaceColor','b');
    % set(h(2),'FaceColor','k');
    
    %Ylabel
    ylabel('\DeltaF/F (%)');
    %ylim,other stuff
    set(gcf,'Color',[1 1 1])
    set(gca,'linewidth',1.6,'FontSize',13,'FontName','Arial')
    
    title({[fpObj.idvData(numMouse).Description ' Before vs After Lick'];...
        ['Range = ' num2str(inspectRange) ' sec']}, 'FontSize',6)
    
    set(gca,'box','off')
    set(gca,'TickDir','out'); % The only other option is 'in'
    
    
    
    %     export_fig([fpObj.idvData(numMouse).Description ' before vs after Lick'],'-eps','-jpg')
    if saveFigures =='y' || saveFigures =='Y'

            print(gcf, '-painters', '-depsc', [fpObj.idvData(numMouse).Description ' barGraphs.svg'])
            print(gcf, '-painters', '-depsc', [fpObj.idvData(numMouse).Description ' barGraphs.pdf'])

    end
end
%% total mean
mOmDff = mean(meanDffStack);
ste_mOmDff = std(meanDffStack)/sqrt(totalMouseNum);
mOmNoneDff = mean(meanNonDffStack);
ste_mOmNoneDff = std(meanNonDffStack)/sqrt(totalMouseNum);
%
figure('visible','on');
hold on
b_1 = bar(1,mOmDff,'FaceColor',[1 1 1],'LineWidth',1.6);
b_2 = bar(2,mOmNoneDff,'FaceColor',[1 1 1],'LineWidth',1.6);
errorbar(1:2,[mOmDff mOmNoneDff],[ste_mOmDff ste_mOmNoneDff],'.','Color','k','linewidth',1.3)
ylim([-1 7]);

% hold on
% h2 = barwitherr(ste_mOmNoneDff, mOmNoneDff)
p = ranksum(meanDffStack,meanNonDffStack)
sigstar([1,2],p)

set(gca,'xtick',[])
Labels = {'Bout', 'NonBout'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
% set(h(1),'FaceColor','b');
% set(h(2),'FaceColor','k');

%Ylabel
ylabel('\DeltaF/F (%)');
%ylim,other stuff
set(gcf,'Color',[1 1 1])
set(gca,'linewidth',1.6,'FontSize',13,'FontName','Arial')
set(gca,'box','off')
set(gca,'TickDir','out'); % The only other option is 'in'



% export_fig([fpObj.idvData(numMouse).Description ' meanOfmean'],'-eps','-jpg')
if saveFigures =='y' || saveFigures =='Y' 
    print(gcf, '-painters', '-depsc', ['boutVSnonBout.svg'])
    print(gcf, '-painters', '-depsc', ['boutVSnonBout.pdf'])
end
%%  mean of mean of mean  

% mOmBeforeLickDff
mOmOmBeforeLickDff= mean(mOmBeforeLickDff);
mOmOmAfterLickDff = mean(mOmAfterLickDff);
%ste of mean
ste_mOmOmBeforeLickDff = std(mOmBeforeLickDff)/sqrt(size(mOmBeforeLickDff,2));
ste_mOmOmAfterLickDff = std(mOmAfterLickDff)/sqrt(size(mOmBeforeLickDff,2));




figure
hold on

b_1 = bar(1,mOmOmBeforeLickDff,'FaceColor',[0 .6 .6],'LineWidth',1.6);
b_2 = bar(2,mOmOmAfterLickDff,'FaceColor',[0 0.2 0.2],'LineWidth',1.6);
errorbar(1:2,[mOmOmBeforeLickDff mOmOmAfterLickDff],[ste_mOmOmBeforeLickDff ste_mOmOmAfterLickDff],'.','Color','k','linewidth',1.3)
p = ranksum(mOmBeforeLickDff',mOmAfterLickDff');
sigstar([1,2],p)

set(gca,'xtick',[])
Labels = {'Before Lick', 'After Lick'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
% set(h(1),'FaceColor','b');
% set(h(2),'FaceColor','k');

%Ylabel
ylabel('\DeltaF/F (%)');
%ylim,other stuff
set(gcf,'Color',[1 1 1])
set(gca,'linewidth',1.6,'FontSize',13,'FontName','Arial')

title({[fpObj.groupInfo{1,1} ' Before vs After Lick: mOmOm'];...
    ['Range = ' num2str(inspectRange) ' sec']}, 'FontSize',10)

set(gca,'box','off')
set(gca,'TickDir','out'); % The only other option is 'in'



% %     export_fig([fpObj.idvData(numMouse).Description ' before vs after Lick'],'-eps','-jpg')
% if saveFigures =='y' || saveFigures =='Y'
%     saveas(gcf,[fpObj.groupInfo{1,1} '  before vs after Lick mOmOm.jpg'])
%     saveas(gcf,[fpObj.groupInfo{1,1} '  before vs after Lick mOmOm.svg'])
% end

%%
% cd ..

end
