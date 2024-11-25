%% Spike simulation
% foodMean = 4.01; foodStd = 0.6*3.32;
waterMean = 0.752; waterStd = 0.0562*3.32;
% Mean = 0.5; Std = 0.1;

timeLimit = 10;
numSamples = 3;
allSpikes = cell(numSamples, 1);

for i = 1:numSamples
    spike = 0; spikes = [];
    % Simulate swallowing pattern for food
    while spike < timeLimit
        % Generate swallowing time
%         interval = normrnd(foodMean, foodStd);
        interval = normrnd(waterMean, waterStd);
%         interval = normrnd(Mean, Std);

        % Update total time
        spike = spike + interval;
        if spike > timeLimit
            break;
        end
        % Store the time
        if size(spikes,2) < 4
            spikes = [spikes spike];
        end
    end
    allSpikes{i} = spikes;
end

examRange = [-15 60];
samplingRate= 102;
time = linspace(examRange(1), examRange(2), (abs(examRange(2))+abs(examRange(1)))*samplingRate+1);

modeledCalcium = [];
for i = 1:numSamples
    optimizedCalcium = calciumModel(fittedParams, allSpikes{i}, time, 10);
    
    %random walk noise
    stdDevIncrements = 0.0121;%std(diff(actualCalcium)); %0.0121 for PB Pdyn citric, 0.0103 for SFO Nos1 citric
    noise = zeros(size(optimizedCalcium));
    % Generate random walk noise
    for i = 2:length(optimizedCalcium)
        noise(i) = noise(i-1) + stdDevIncrements * randn;
    end
    
    optimizedCalcium = optimizedCalcium + noise;
    modeledCalcium = [modeledCalcium; optimizedCalcium];
end

a=[a; mean(modeledCalcium)];

%%
figure;
% plot(time, mean(modeledCalcium), 'r-', 'LineWidth', 2);
mseb(time,mean(modeledCalcium),std(modeledCalcium)./sqrt(size(modeledCalcium,1)));
hold on;
% legend('Actual Calcium', 'Fitted Model');
% ylim([-3 1])
xlabel('Time');
ylabel('Calcium Level');
set(gca,'linewidth',1.6,'FontSize',15,'FontName','Arial')
set(gca, 'box', 'off')
set(gcf,'Color',[1 1 1])
set(gca,'TickDir','out'); % The only other option is 'in'

hold off
%%
figure;

for numMouse =1:size(allSpikes,1)
    spikes = allSpikes{numMouse};
    
    %Shading licking as dark red
    for i = 1:size(spikes,2)
        line([spikes(i) spikes(i)], [numMouse-1 numMouse],'Color','r')
    end
%     xlim([0 1800]);
    
%     %Shading bout as light red
%     for i = 1:size(boutIdx,1)
%         r = patch([boutIdx(i,1) boutIdx(i,2) boutIdx(i,2) boutIdx(i,1)], [0 0 1 1],...
%             [1,0,0]);
%         set(r, 'FaceAlpha', 0.2,'LineStyle','none');
%         uistack(r,'up')
%     end
end

xlim([0 10])