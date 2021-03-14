%% Calculate NRMS for fits
clear all;
load('./Data/exportDataFits.mat');

for ii = 1:length(exportData)
    x = exportData{ii}.detuning;
    y = exportData{ii}.Ts1550;
    y_fit = exportData{ii}.Ts1550fit;

    [NRMSE(ii)] = calculateNRMSE(x, y, y_fit');
end

%%
load('./Data/exportDataModel.mat');

for ii = 1:length(exportData)
    x = exportData{ii}.detuning;
    y = exportData{ii}.Ts1550;
    y_fit = exportData{ii}.Ts1550fit;

    [NRMSE_model(ii)] = calculateNRMSE(x, y, y_fit');
end

%% Plots

myFig = figure(); 
plot([1:3, 5:7, 9], NRMSE([1:3, 5:7, 9])*100, '.'); hold on;
plot(NRMSE_model*100, '.');

xlabel('Dataset #');
ylabel('NRMSE (%)');
xlim([0 10]);
xticks([1:2:9]);

saveas(myFig, './final_plots/NRMSE.fig');
