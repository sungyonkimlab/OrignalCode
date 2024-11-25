clear all; close all

%%
examRange = [-15 60];
samplingRate = 102;

c = {processedData{8,2}}; %water
c = {c{1,1} processedData{15,2}}; %mechano

d = {dFFStack_inRange(:,8)}; % water
d = {d{1,1} dFFStack_inRange(:,15)}; % meechano

%% PB Pdyn
spikeTrains = c;
calciumData = d;
time = linspace(examRange(1), examRange(2), size(d{1},1));
timeVectors = {time, time};
fMax = {}; calciumDatas = {}; actualCalcium = {};
for idx = 1:size(d,2)
    actualCalcium{idx} = d{idx}'-mean(d{idx}(1:abs(examRange(1))*samplingRate));
    fMax{idx} = max(actualCalcium{idx});
    calciumDatas{idx} = actualCalcium{idx};
end
fMax = max(fMax{:});

%% SFO thirst
spikeTrains = c;
calciumData = d; 
time = linspace(examRange(1), examRange(2), size(d{1},1));
timeVectors = {time, time};
fMax = {}; 
calciumDatas = {}; actualCalcium = {};
for idx = 1:size(d,2)
    actualCalcium{idx} = -(d{idx}'-mean(d{idx}(1:abs(examRange(1))*samplingRate)));
    fMax{idx} = max(smooth(actualCalcium{idx},samplingRate));
    fMax{idx} = max(actualCalcium{idx});
    calciumDatas{idx} = actualCalcium{idx};
end
fMax = max(fMax{:});

%%
% initial parameter estimate
initialParams = [1, 0.1, 10, 1, 1, 0, 0]; % [tauRise, tauDecay1, tauDecay2, r, k, C_half, delay]
lb = [0, 0, 0, 0,  0, -Inf, 0]; % lower bound
ub = [Inf, Inf, Inf, Inf, Inf, Inf, Inf]; % upper bound

% lsqnonlin options
options = optimoptions('lsqnonlin', 'Display', 'iter','MaxFunctionEvaluations', 1600, 'FunctionTolerance', 1e-6);

% lsqnonlin
fittedParams = lsqnonlin(@(p)calcErrorAcrossDatasets(p, spikeTrains, calciumDatas, timeVectors, fMax), initialParams, lb, ub, options);

% diplay tuned parameters [tauRise, tauDecay1, tauDecay2, r, A, k, C_half, l]
disp('Fitted Parameters:');
disp(['tauRise: ' num2str(fittedParams(1))]);
disp(['tauDecay1: ' num2str(fittedParams(2))]);
disp(['tauDecay2: ' num2str(fittedParams(3))]);
disp(['r: ' num2str(fittedParams(4))]);
disp(['k: ' num2str(fittedParams(5))]);
disp(['C_half: ' num2str(fittedParams(6))]);
disp(['delay: ' num2str(fittedParams(7))]);
% disp(['fMax: ' num2str(fittedParams(8))]);

%%
optimizedCalcium = {};
% optimizedC = {};
% optimizedF = {};
for idx = 1:size(d,2)
    optimizedCalcium{idx} = calciumModel(fittedParams, spikeTrains{idx}, timeVectors{idx}, fMax);
end

%%
% result visualization
figure;
plot(time, -actualCalcium{2}, 'b-', 'LineWidth', 2);
hold on;
plot(time, -optimizedCalcium{2}, 'r-', 'LineWidth', 2);
% legend('Actual Calcium', 'Fitted Model');
title('Calcium Data and Fitted Model');
xlabel('Time');
ylabel('Calcium Level');

%%
t=0:0.01:10;
visKernel = (fittedParams(4)*exp(-t/fittedParams(2)) + exp(-t/fittedParams(3))) .* (1-exp(-t/fittedParams(1)));
figure;
plot(t,visKernel)

%%
x=0:0.01:ceil(5*max(visKernel));
sigmoid = 1 ./ (1 + exp(-fittedParams(5) * (x - fittedParams(6))));
xSwallow = 0:1*max(visKernel):5*(max(visKernel));
% xSwallow = [0];
% for i = 0:4
%     xSwallow = [xSwallow max(cModel(fittedParams,0:i,timeVectors{1}))];
% end
ySwallow = 1 ./ (1 + exp(-fittedParams(5) * (xSwallow - fittedParams(6))));

figure;
plot(x,sigmoid)
hold on
scatter(xSwallow, ySwallow,12)

xlim([0, max(xSwallow)])
set(gca,'linewidth',1.6,'FontSize',15,'FontName','Arial')
set(gca, 'box', 'off')
set(gcf,'Color',[1 1 1])
set(gca,'TickDir','out'); % The only other option is 'in'
yticks(0:0.2:1)
ylim([0 1])

%% citric acid validataion
% dFFStack_inRange
k=20;
spikes = processedData{k,2};
actualCalcium = -(dFFStack_inRange(:,k)'-mean(dFFStack_inRange(1:abs(examRange(1))*samplingRate,k)));
% fMax = max(smooth(actualCalcium,samplingRate));

examRange = [-15 60];
time = linspace(examRange(1), examRange(2), size(d{1},1));
fMax = max(max(d{:}));

optimizedCalcium = calciumModel(fittedParams, spikes, time, fMax);

figure;
plot(time, -actualCalcium, 'b-', 'LineWidth', 2);
hold on;
plot(time, -optimizedCalcium, 'r-', 'LineWidth', 2);
hold on
% title('Calcium Data and Fitted Model');
xlabel('Time');
ylabel('Calcium Level');
set(gca,'linewidth',1.6,'FontSize',15,'FontName','Arial')
set(gca, 'box', 'off')
set(gcf,'Color',[1 1 1])
set(gca,'TickDir','out'); % The only other option is 'in'

explainedVar = 1- var(actualCalcium-optimizedCalcium)/var(actualCalcium)