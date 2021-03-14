%% Find all data paths

clear all; close all;

addpath('./helpers');

color_palette;

aux_plots               = 0;
data_stamp_target       = '004';
data_stemp_target_num   = str2double(data_stamp_target);

dinfo = dir('./ExperimentalData');
dinfo(ismember( {dinfo.name}, {'.', '..'})) = []; 

for ii = 1:length(dinfo)
    data_stamp(ii) = getparam(dinfo(ii).name, 'DataPowerAndT_');
end

directory   = [dinfo( data_stamp == data_stemp_target_num ).folder '\' dinfo( data_stamp == data_stemp_target_num ).name];

dinfo = dir(directory);
dinfo(ismember( {dinfo.name}, {'.', '..', 'Processed'})) = []; 

for ii = 1:length(dinfo)
   temperature(ii) = getparam(dinfo(ii).name, '_T_'); 
end

[temperature, ind] = sort(temperature);
dinfo = dinfo(ind);

%%
tt = 1;
directory   = [ dinfo(tt).folder '/' dinfo(tt).name ];
saveDR      = [ './final_plots/' ];
mkdir(saveDR);

%% Load data

dd = filefun([ directory '/power_sweeps/LNOI15_*.mat']);

for ii = 1:length(dd)
    rep_no(ii) =  getparam(dd{ii}, 'rep_');
end

[rep_no, ind] = sort(rep_no);
dd = dd(ind);

%% Load experimental parameters for further processing

repetition_no = 1;
gg = repetition_no;

load(dd{1});
power_points = length(data);

load(dd{gg});

for ii = 1:power_points

   power1550(ii)                                = data{ii}.input_PM_1550 * cal{ii}.W_to_uW_1550_input_in_fiber;
   attens(ii)                                   = cal{ii}.input_atten;
   data_local(ii).wls                           = data{ii}.wls;
   
   SH_index                 = find( strcmp(cal{ii}.daq_info, 'SHG') );
   Fundamental_index        = find( strcmp(cal{ii}.daq_info, 'Fundamental') );
   
   data_local(ii).Vs1550                        = data{ii}.Vs(:,Fundamental_index);
   data_local(ii).Vs780                         = data{ii}.Vs(:,SH_index);
   data_local(ii).backgroundV                   = cal{ii}.backgroundV;
   data_local(ii).backgroundV_780               = cal{ii}.backgroundV_780;
   data_local(ii).output_atten1550              = cal{ii}.output_atten1550;
   data_local(ii).ND_filter                     = cal{ii}.ND_filter;

   data_local(ii).V_to_uW_1550_output_in_fiber  = cal{ii}.V_to_uW_1550_output_in_fiber;
   data_local(ii).V_to_uW_780_output_in_fiber   = cal{ii}.V_to_uW_780_output_in_fiber;

   
   data_local(ii).Power780_in_fiber             = (data_local(ii).Vs780 - cal{ii}.backgroundV_780) * data_local(ii).V_to_uW_780_output_in_fiber;
   data_local(ii).Power1550_in_fiber            = (data_local(ii).Vs1550 - cal{ii}.backgroundV) * cal{ii}.V_to_uW_1550_output_in_fiber;
   
   data_local(ii).fiber_to_chip_1550            = sqrt( abs(mean(data_local(ii).Power1550_in_fiber)) / power1550(ii) );
%    data_local(ii).fiber_to_chip_1550            = coupling_efficiency;
   data_local(ii).Power1550_on_chip             = (data_local(ii).Vs1550 - cal{ii}.backgroundV) * cal{ii}.V_to_uW_1550_output_in_fiber / data_local(ii).fiber_to_chip_1550;

   power1550_on_chip(ii)                        = power1550(ii)* data_local(ii).fiber_to_chip_1550;
   
end

%% Overwrite coupling efficiency with average value
calculateCouplingEfficiency;
for ii = 1:power_points
	data_local(ii).fiber_to_chip_1550            = coupling_efficiency;
    data_local(ii).Power1550_on_chip             = (data_local(ii).Vs1550 - cal{ii}.backgroundV) * cal{ii}.V_to_uW_1550_output_in_fiber / data_local(ii).fiber_to_chip_1550;

    power1550_on_chip(ii)                        = power1550(ii)* data_local(ii).fiber_to_chip_1550;
end

%% Sort data wrt. power
[attens, ind]               = sort(attens, 'descend');
power1550                   = power1550(ind);
data_local                  = data_local(ind);
power1550_on_chip           = power1550_on_chip(ind);

%% Extract data and fits for class processing

load('./Data/lineshape_plot_data.mat');

indices = 29:1:37;
myFig = figure();
ranges      = linspace(0.8, 1.1750, length(indices)+1 )*1e10;
iteration = 1;

for ii = indices
    
    step    = (ranges(1) + ranges(iteration + 1) )/2 + 0.1e10;
    
    range   = ranges(iteration);
    filter  = abs(plot_data{ii}.detuning)<range/2;
    
    detuning_local  = plot_data{ii}.detuning(filter);
    Ts1550_loc      = plot_data{ii}.Ts1550(filter);
    SHG_loc         = plot_data{ii}.SHG(filter);
    
    % Get background
    toFitX  = [mean(detuning_local(1:10)), mean(detuning_local(end-10:end))];
    toFitY  = [mean(Ts1550_loc(1:10)), mean(Ts1550_loc(end-10:end))];
    FF      = polyfit(toFitX, toFitY, 1); 
    bgfit   = FF(1) * detuning_local + FF(2);
    
    Ts1550_loc = Ts1550_loc ./bgfit;
    
    if ~isempty(plot_data{ii}.Ts1550_fitted) 
        Ts1550_fit_loc  = plot_data{ii}.Ts1550_fitted(filter)' ./ bgfit';
        SHG_fit_loc     = plot_data{ii}.SHG_fitted(filter);
    end
    
    plot(detuning_local + (iteration-1)*step, Ts1550_loc, '.', 'MarkerSize', datapoint_size, 'Color', [hex2rgb(color.faint_orange), 1]); hold on;
    if ~isempty(plot_data{ii}.Ts1550_fitted)
        plot(detuning_local + (iteration-1)*step, Ts1550_fit_loc, '-', 'Color', [hex2rgb(color.orange), 1], 'LineWidth', line_width2); hold on;
    end

    plot(detuning_local + (iteration-1)*step, SHG_loc, '.', 'MarkerSize', datapoint_size, 'Color', [hex2rgb(color.faint_blue), 1]); hold on;
    if ~isempty(plot_data{ii}.SHG_fitted)
        plot(detuning_local + (iteration-1)*step, SHG_fit_loc, '-', 'Color', [hex2rgb(color.blue), 1], 'LineWidth', line_width2); hold on;
    end

    exportData{iteration}.detuning  = detuning_local + (iteration-1)*step;
    exportData{iteration}.Ts1550    = Ts1550_loc;
    exportData{iteration}.Ts1550fit = Ts1550_fit_loc;

    iteration = iteration + 1;

end

save('./Data/exportDataModel.mat', 'exportData');
saveas(myFig, './final_plots/allFitsModel.fig');

%% Extract data and fits for class processing

load('./Data/lineshape_plot_data_fits.mat');

indices = 29:1:37;
myFig = figure();
ranges      = linspace(0.8, 1.1750, length(indices)+1 )*1e10;
iteration = 1;

for ii = indices
    
    step    = (ranges(1) + ranges(iteration + 1) )/2 + 0.1e10;
    
    range   = ranges(iteration);
    filter  = abs(plot_data{ii}.detuning)<range/2;
    
    detuning_local  = plot_data{ii}.detuning(filter);
    Ts1550_loc      = plot_data{ii}.Ts1550(filter);
    SHG_loc         = plot_data{ii}.SHG(filter);
    
    % Get background
    toFitX  = [mean(detuning_local(1:10)), mean(detuning_local(end-10:end))];
    toFitY  = [mean(Ts1550_loc(1:10)), mean(Ts1550_loc(end-10:end))];
    FF      = polyfit(toFitX, toFitY, 1); 
    bgfit   = FF(1) * detuning_local + FF(2);
    
    Ts1550_loc = Ts1550_loc ./bgfit;
    
    if ~isempty(plot_data{ii}.Ts1550_fitted) 
        Ts1550_fit_loc  = plot_data{ii}.Ts1550_fitted(filter)' ./ bgfit';
        SHG_fit_loc     = plot_data{ii}.SHG_fitted(filter);
    end
    
    plot(detuning_local + (iteration-1)*step, Ts1550_loc, '.', 'MarkerSize', datapoint_size, 'Color', [hex2rgb(color.faint_orange), 1]); hold on;
    if ~isempty(plot_data{ii}.Ts1550_fitted)
        plot(detuning_local + (iteration-1)*step, Ts1550_fit_loc, '-', 'Color', [hex2rgb(color.orange), 1], 'LineWidth', line_width2); hold on;
    end

    plot(detuning_local + (iteration-1)*step, SHG_loc, '.', 'MarkerSize', datapoint_size, 'Color', [hex2rgb(color.faint_blue), 1]); hold on;
    if ~isempty(plot_data{ii}.SHG_fitted)
        plot(detuning_local + (iteration-1)*step, SHG_fit_loc, '-', 'Color', [hex2rgb(color.blue), 1], 'LineWidth', line_width2); hold on;
    end

    exportData{iteration}.detuning  = detuning_local + (iteration-1)*step;
    exportData{iteration}.Ts1550    = Ts1550_loc;
    exportData{iteration}.Ts1550fit = Ts1550_fit_loc;

    iteration = iteration + 1;

end

save('./Data/exportDataFits.mat', 'exportData');
saveas(myFig, './final_plots/allFitsFits.fig');

