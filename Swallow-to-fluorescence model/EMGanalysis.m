clear all; close all;

%% Read data
rawData = TDTbin2mat('D:\MY lab\Projects\_rapid thirst satiation\Fiber photometry\SFO Nos1 GCaMP 7th\2024-03-12~27 SFO Nos1 GCaMP swallow, probing after pancuronium\_EMG\240312_Panc_Nos1_4345');%file location

%% Construct data structure for analysis
channelName = {"EMG1"};
channelData = double(rawData.streams.("EMG1").data)';
channelFs = rawData.streams.("EMG1").fs;
channelTime = ((0:length(channelData)-1) / channelFs)';

data = struct('channelName', channelName, 'channelData', channelData, 'channelFs', channelFs, 'channelTime', channelTime);
eventTime = rawData.epocs.Ep1_.onset;

%% Rectify and filter the signal (EMG)
order = 2;
fcLow = 5;
fcHigh = 0.05;

[b, a]  = butter(order, fcLow/(channelFs/2), 'low');
[bb, aa] = butter(order, fcHigh/(channelFs/2), 'high');
rectifiedSignal = abs(channelData);

filteredSignal = filtfilt(b,a, rectifiedSignal);
filteredSignal = filtfilt(bb,aa,filteredSignal);

channelData = filteredSignal;

%% Plot original and fltered traces
t = (0:length(channelData)-1)/channelFs;
legend('Input Data','Filtered Data')
yyaxis left;
plot(t, channelData, 'b-'); % First channel in blue
ylabel("Original trace");
yyaxis right;
plot(t, filteredSignal, 'r-'); % Second channel in red
ylabel("Filtered trace");
xlabel('Time');

%% Detect onsets around events
%set parameters for onset detection
examRange = [-15 60];
% baselinePeriod = [-10,-2];
thresholdMultiplier = 5; % set threshold (x baseline sd)
minInterval = 0.3; % minimal interval between onsets

processedData = cell(size(eventTime, 1), 5);
for eventNum = 1:size(eventTime, 1)
    eventTimeIdx = round(eventTime(eventNum) * channelFs);
    windowStartIdx = eventTimeIdx + round(examRange(1) * channelFs);
    windowEndIdx = eventTimeIdx + round(examRange(2) * channelFs);
    sampledData = channelData(windowStartIdx:windowEndIdx);
    
    % Calculate threshold
    baselineStartIdx = windowStartIdx;
    baselineEndIdx = windowEndIdx;
%     baselineStartIdx = round((eventTime(eventNum) + baselinePeriod(1)) * channelFs);
%     baselineEndIdx = round((eventTime(eventNum) + baselinePeriod(2)) * channelFs);
    baselineData = channelData(baselineStartIdx:baselineEndIdx);
    meanBaseline = mean(baselineData);
    stdBaseline = std(baselineData);
    threshold = meanBaseline + thresholdMultiplier * stdBaseline;
    
    % Detect onsets
    onsetTimes = [];
    lastOnsetIdx = 0; % Initialize to start of detection window
    belowThreshold = true; % Initialize the below-threshold flag
    
    for i = 1:size(sampledData,1)
        if belowThreshold && sampledData(i) > threshold && (i - lastOnsetIdx) > round(minInterval * channelFs)
            onsetTime = (windowStartIdx + i - 1) / channelFs; % Convert index to time
            onsetTimes(end + 1) = onsetTime;
            lastOnsetIdx = i;
            belowThreshold = false;
        elseif sampledData(i) < threshold
            belowThreshold = true; % Set the flag when signal goes below threshold
        end
    end
    numOnsets = length(onsetTimes);
    
    % Store the results in the corresponding row of the cell array
    processedData{eventNum, 1} = threshold;
    if ~isempty(onsetTimes)
        processedData{eventNum, 2} = onsetTimes - eventTime(eventNum); % latency
    end
    processedData{eventNum, 3} = sampledData;
    processedData{eventNum, 4} = onsetTimes;
    processedData{eventNum, 5} = numOnsets;
end

%% Plot EMG trace and mark eventTime, baseline, onsets, and threshold
for eventNum = 1:size(processedData, 1)
    sampledData = processedData{eventNum, 3};
    onsetTime = processedData{eventNum, 2};
    threshold = processedData{eventNum, 1};
    
    figure;
    plot(linspace(examRange(1),examRange(2),size(sampledData,1)), sampledData);
    hold on;
    ylabel('EMG amplitude');
    xlabel('Time (seconds)');
    
    yline(threshold, 'k--', 'LineWidth', 1.5);
    
    for i = 1:length(onsetTime)
%         xline(onsetTime(i), 'r-', 'LineWidth', 1.5);
%         text(onsetTime(i), max(sampledData), sprintf(' # %d', i), 'Color', 'red', 'VerticalAlignment', 'bottom', 'FontSize',10);
    end
    
    set(gca,'linewidth',1.6,'FontSize',15,'FontName','Arial')
    set(gca, 'box', 'off')
    set(gca,'TickDir','out'); % The only other option is 'in'
    
    xlim(examRange)
    
    hold off;
    saveas(gcf,[num2str(eventNum) '.jpeg'])
end