%% current consumption measured with KEYSIGNT N6705B

clc
clear all
filename = ['partial_analysis/dlog10.csv']; % 10, 11

M = csvread(filename,7,0);
sampleTime = 0.00008192;  % 8.192e-005
data.time = M(:,1).*sampleTime;   
data.voltage = M(:,2);
data.current = M(:,3) * 1000;
data.sampleTime = sampleTime;
clear M  filename sampleTime


%% dataNode2.time,smooth(dataNode2.current,5000),dataNode3.time,smooth(dataNode3.current,0.005,'rloess')

data.current_smoothed = smooth(data.current, 150);
data.power = data.voltage .* data.current;
data.power_smoothed = smooth(data.power, 100);

figure
plot(data.time, data.current, 'Color',[0 0 0]+0.8);
hold on
plot(data.time, data.current_smoothed);
xlabel('Time (s)')
ylabel('Current (mA)')
legend('Current', 'Average', 'Location','best'); 



figure
plot(data.time, data.power, 'Color',[0 0 0]+0.8);
hold on
grid on
plot(data.time, data.power_smoothed, 'linewidth', 1.5);
xlabel('Time (s)')
ylabel('Power (mW)')
% ylim([0, 3.5])
legend('Power', 'Average', 'Location','best'); 


%%

i_start = 15492;
i_end = 23194; % size(data.current,1);

plot_time = data.time(1:i_end-i_start+1);
plot_power = data.power(i_start:i_end);
plot_power_smoothed = data.power_smoothed(i_start:i_end);

figure
plot(plot_time, plot_power, 'Color',[0 0 0]+0.8);
hold on
grid on
plot(plot_time, plot_power_smoothed, 'linewidth', 1.5);
xlabel('Time (s)')
ylabel('Power (mW)')
% xlim([0, 2])
legend('Power', 'Average', 'Location','best'); 
set(gca, 'fontsize', 14)
figurepalette