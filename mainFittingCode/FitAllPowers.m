%%
clear all; close all;
addpath('./helpers');

modeParameters;

%% Load data

load('./Data/LNOI15_dev36_rep_1.mat');

fiber_to_chip780 = 0.11;
fiber_to_chip1550 = 0.26;

%%%%%%%%%%%%%%%%%%%%%%
%% FITTING SECTION %%%
%%%%%%%%%%%%%%%%%%%%%%


datestamp = datestr(now,30);
mkdir(['./FitLineshapes' datestamp '/']);

for kk = 38:-1:14

data_loc  = data{kk};
cal_loc   = cal{kk};

%

Ts  = (data_loc.Vs(:,1) - cal_loc.backgroundV_780) * cal_loc.V_to_uW_780_output_in_fiber / fiber_to_chip780 * 1e-6;
wls = data_loc.wls;
P_pump = data_loc.input_PM_1550 * cal_loc.W_to_uW_1550_input_in_fiber * fiber_to_chip1550 * 1e-6;
Ts1550 = (data_loc.Vs(:,2) - cal_loc.backgroundV) * cal_loc.V_to_uW_1550_output_in_fiber / fiber_to_chip1550 * 1e-6; 
TsMZI  = data_loc.Vs(:,3);

% wl_filter = wls > 1549.407 & wls < 1549.42;
% wl_filter = wls > 1549.394 & wls < 1549.431;
wl_filter = wls > 1549.405 & wls < 1549.44;    % For power 35

Ts = Ts(wl_filter);
wls = wls(wl_filter);
Ts1550 = Ts1550(wl_filter)';
TsMZI = TsMZI(wl_filter);

% MZI

MZI1550_interp   = interp1(1:length(wls), TsMZI, 1:0.1:length(wls));
wls_interp      = interp1(1:length(wls), wls, 1:0.1:length(wls));

[wls_out Ts780_norm_out delta_lambda_calibration] = mzcalibration_tmp( wls_interp, MZI1550_interp, MZI1550_interp', 0.0, 'MARKO');

wls         = interp1(1:length(wls_out), wls_out, 1:10:length(wls_out));

% % Narrow filtering
% wl_filter = wls > 1549.409 & wls < 1549.418;
% % wl_filter = wls > 1549.394 & wls < 1549.431;

% Find center wavelength
Ts1550 = Ts1550/max(Ts1550);

[sh_max, max_ind] = max(smooth(1-Ts1550));
center_wl = mean(wls((1-Ts1550)/max(1-Ts1550) > 0.7));

% Narrow filtering
wl_window = 0.02;
wl_filter = (wls > (center_wl-wl_window/2)) & (wls < (center_wl+wl_window/2));

Ts = Ts(wl_filter);
wls = wls(wl_filter);
Ts1550 = Ts1550(wl_filter);


wls_to_fit = [mean(wls(1:10)) mean(wls(end-10:end))];
Ts1550_to_fit = [mean(Ts1550(1:10)) mean(Ts1550(end-10:end))];

bgFitP = polyfit(wls_to_fit, Ts1550_to_fit, 1);
bgFit = wls.*bgFitP(1) + bgFitP(2);

figure(); 
subplot(211);
plot(wls, Ts1550); hold on;
plot(wls, bgFit);

% Correct for bg
Ts1550 = Ts1550 ./ bgFit;
subplot(212);
plot(wls, Ts1550);

% Find center wavelength
[sh_max, max_ind] = max(smooth(Ts));
center_wl = mean(wls(Ts/max(Ts) > 0.5))*1e-9;

Ts = Ts(1:3:end);
wls = wls(1:3:end)*1e-9;
Ts1550 = Ts1550(1:3:end);

% lambda_c = 1549.41385*1e-9;     % For power
% lambda_c = 1549.41385*1e-9;

% Try fitting
clear FitParameters;
F = @(FitParameters, x)SHG_amplitude(FitParameters, x);

wavelength_a     = center_wl; %lambda_c; % 1.549414e-06; % lambda_c;
Qa               = Q_a; %0.9e6; % Q_a;
Qb               = Q_b;
g                = g_experiment/2/pi;
coupling_ratio_a = kappa_a_e/kappa_a;       
coupling_ratio_b = kappa_b_e/kappa_b;       % doesn't impact transmission lieshape

lines_detuning   = 2.691942078220962e+08;
% P_pump           = 2e-4;

FitParameters = [wavelength_a Qa Qb  ...
                g coupling_ratio_a coupling_ratio_b ...
                lines_detuning P_pump];
            
% Fit transmission
F = @(FitParameters, x)SHG_pump_transmission(FitParameters, x);

% what_to_fit = [true false true false false true false true];
what_to_fit = [true false true false false true false true];

opts = statset('nlinfit');
opts.MaxIter = 200;

tic;
NewParameters = nlinfitsome(what_to_fit, wls, Ts1550, F , FitParameters, opts);
toc

% Refine g with the SHG lineshape
F = @(FitParameters, x)SHG_amplitudeNorm(FitParameters, x);

% what_to_fit = [true false true false false true true true];
what_to_fit = [true false true false false true true true];

Ts = Ts' / max(smooth(Ts,3));

tic;
NewParameters = nlinfitsome(what_to_fit, wls, Ts, F , NewParameters, opts);
toc

% Refine fit with transmission
F = @(FitParameters, x)SHG_pump_transmission(FitParameters, x);

% what_to_fit = [true false true true false true true true];
what_to_fit = [true false false true false true false true];

tic;
NewParameters = nlinfitsome(what_to_fit, wls, Ts1550, F , NewParameters, opts);
toc

% Plot fitted data
FitParameters(5) = 0.35;
y_fitted = @(x) SHG_pump_transmission(NewParameters, x);

wls2 = wls;
y2 = y_fitted(wls);

y_starting = @(x) SHG_pump_transmission(FitParameters, x);

wls3 = wls;
y3 = y_starting(wls3);

y_fitted2 = @(x) SHG_amplitudeNorm(NewParameters, x);

wls2 = wls;
y4 = y_fitted2(wls);

y_starting2 = @(x) SHG_amplitudeNorm(FitParameters, x);

wls3 = wls;
y5 = y_starting2(wls3);

fit_fig(kk) = figure();
subplot(121);
plot( wls, Ts1550, 'o' ); hold on;
plot( wls3, y3, '-', 'color', 'green');
plot( wls2, y2, '-', 'color', 'red'); 
grid;
subplot(122);
plot( wls, Ts/max(smooth(Ts)), 'o' ); hold on;
plot( wls3, y5, '-', 'color', 'green');
plot( wls2, y4, '-', 'color', 'red'); 
grid;


subplot(121);
ax = gca;

Qa_fitted = NewParameters(2);
g_fitted = NewParameters(4);
coupling_ratio_a_fitted = NewParameters(5);

text(ax.XLim(1), 0.7, ['Start parameters:' ]);
text(ax.XLim(1), 0.65, ['Qa: ' num2str(FitParameters(2)) ]);
text(ax.XLim(1), 0.6, ['\kappa_e/\kappa: ' num2str(FitParameters(5)) ]);
text(ax.XLim(1), 0.55, ['g: ' num2str(FitParameters(4)) ]);

text(ax.XLim(1), 0.45, ['Fitted:' ]);
text(ax.XLim(1), 0.4, ['Qa: ' num2str(Qa_fitted) ]);
text(ax.XLim(1), 0.35, ['\kappa_e/\kappa: ' num2str(coupling_ratio_a_fitted) ]);
text(ax.XLim(1), 0.3, ['g: ' num2str(g_fitted) ]);
text(ax.XLim(1), 0.25, ['Qae: ' num2str(Qa_fitted / coupling_ratio_a_fitted) ]);

legend('Data', 'Start parameters', 'Fitted');

allFitParameters(kk).NewParameters = NewParameters;

saveas(fit_fig(kk), ['./FitLineshapes' datestamp '/data_' num2str(kk) '.fig']);
saveas(fit_fig(kk), ['./FitLineshapes' datestamp '/data_' num2str(kk) '.png']);

end

%% Save fit parameters

save(['./FitLineshapes' datestamp '/fitData.mat'], 'allFitParameters');

%% Plot fit results

clear g Qa wa ka kae kai Qae Qai Detuning Power_in Qb;
iteration = 1;

for ii = 20:length(allFitParameters)
    
    data_loc  = data{ii};
    cal_loc   = cal{ii};
    
   g(iteration)    = allFitParameters(ii).NewParameters(4);
   Qa(iteration)   = allFitParameters(ii).NewParameters(2);
   Qb(iteration)   = allFitParameters(ii).NewParameters(3);
   
   wa(iteration)   = 2*pi*c/allFitParameters(ii).NewParameters(1);
   ka(iteration)   = wa(iteration)/Qa(iteration);
   
   kae(iteration)  = allFitParameters(ii).NewParameters(5)*ka(iteration);
   kai(iteration)  = ka(iteration) - kae(iteration);
   
   Qae(iteration)  = wa(iteration)/kae(iteration);
   Qai(iteration)  = wa(iteration)/kai(iteration);
   
   Detuning(iteration) = allFitParameters(ii).NewParameters(7);
   Power_in(iteration) = data_loc.input_PM_1550 * cal_loc.W_to_uW_1550_input_in_fiber * fiber_to_chip1550 * 1e-6;
    
   iteration = iteration + 1;
   
end

save(['./FitLineshapes' datestamp '/fitData.mat'], 'allFitParameters', 'g', 'Qa', 'Qb', 'wa', 'ka', 'kae', 'kai', 'Qae', 'Qai', 'Detuning', ...
    'Power_in');

%%

datestamp = '20210313T191602';

cutoff = 9;
cutoff_end = 18;

mean_Qa     = mean(Qa(cutoff:cutoff_end)); 
mean_Qae    = mean(Qae(cutoff:cutoff_end));
mean_Qb     = mean(Qb(cutoff:cutoff_end)); 
mean_g      = mean(g(cutoff:cutoff_end));

all_Qs_fig = figure();
subplot(311);
loglog(Power_in, Qa, 'o-');  hold on;
loglog(Power_in, Qb, 'o-'); legend(['Qa, mean = ' num2str(mean_Qa) ', std.d. = ' num2str(std(Qa(cutoff:end))) ], ...
    ['Qb, mean = ' num2str(mean_Qb) ', std.d. = ' num2str(std(Qb(cutoff:end))) ]);
subplot(312);
loglog(Power_in, Qae, 'o-'); legend(['Qae, mean = ' num2str(mean_Qae) ', std.d. = ' num2str(std(Qae(cutoff:end))) ] );
subplot(313); 
loglog(Power_in, g, 'o-'); legend(['g, mean = ' num2str(mean_g) ', std.d. = ' num2str(std(g(cutoff:end))) ] );
xlabel('P_{in} on chip');

saveas(all_Qs_fig, ['./FitLineshapes' datestamp '/AllQs_gs.fig']);
saveas(all_Qs_fig, ['./FitLineshapes' datestamp '/AllQs_gs.png']);

%% Plot detuning
mean_k = mean(ka);

allDetunings_fig = figure();
semilogx(Power_in, Detuning, 'o-'); 
yline(mean_k);
yline(-mean_k);
grid; ylim([-2e9 2e9]);
legend(['Detuning, mean = ' num2str(mean(Detuning(18:end))/mean_k) '\kappa, std.d. = ' num2str(std(Detuning(18:end))/mean_k) '\kappa' ], ...
    '+\kappa_a', '-\kappa_a');

saveas(allDetunings_fig, ['./FitLineshapes' datestamp '/AllDetunings.fig']);
saveas(allDetunings_fig, ['./FitLineshapes' datestamp '/AllDetunings.png']);

%% Plot theory lineshapes on top of the experimental data

all_lineshapes = figure();
iteration = 1;

modeParameters;

for kk = 1:1:length(data)

% Load single data

data_loc  = data{kk};
cal_loc   = cal{kk};

%
fiber_to_chip780 = 0.09;
fiber_to_chip1550 = 0.26;

Ts  = (data_loc.Vs(:,1) - cal_loc.backgroundV_780) * cal_loc.V_to_uW_780_output_in_fiber / fiber_to_chip780 * 1e-6;
wls = data_loc.wls*1e-9;
P_pump = data_loc.input_PM_1550 * cal_loc.W_to_uW_1550_input_in_fiber * fiber_to_chip1550 * 1e-6;
Ts1550 = (data_loc.Vs(:,2) - cal_loc.backgroundV) * cal_loc.V_to_uW_1550_output_in_fiber / fiber_to_chip1550 * 1e-6; 
TsMZI  = data_loc.Vs(:,3);

wl_filter = wls > 1549.35*1e-9 & wls < 1549.47*1e-9;    % For power 35

Ts = Ts(wl_filter);
wls = wls(wl_filter);
Ts1550 = Ts1550(wl_filter)';
TsMZI = TsMZI(wl_filter);

% MZI

MZI1550_interp   = interp1(1:length(wls), TsMZI, 1:0.1:length(wls));
wls_interp      = interp1(1:length(wls), wls, 1:0.1:length(wls));

[wls_out Ts780_norm_out delta_lambda_calibration] = mzcalibration_tmp( wls_interp*1e9, MZI1550_interp, MZI1550_interp', 0.0, 'MARKO');

wls         = interp1(1:length(wls_out), wls_out, 1:10:length(wls_out))*1e-9;
Ts1550      = Ts1550/max(Ts1550);

% Find center wavelength
[sh_max, max_ind] = max(smooth(1-Ts1550));
center_wl = mean(wls((1-Ts1550)/max(1-Ts1550) > 0.7));

% Narrow filtering
wl_window = 0.02*1e-9;
if kk >= 28 && kk <=38
    wl_filter = (wls > (center_wl-wl_window/2-0.005*1e-9)) & (wls < (center_wl+wl_window/2+0.000*1e-9));
elseif kk > 38
    wl_filter = (wls > (center_wl-wl_window/2-0.01*1e-9)) & (wls < (center_wl+wl_window/2+0.01*1e-9));
else
    wl_filter = (wls > (center_wl-wl_window/2-0.00*1e-9)) & (wls < (center_wl+wl_window/2+0.000*1e-9));
end

Ts = Ts(wl_filter);
wls = wls(wl_filter);
Ts1550 = Ts1550(wl_filter);

wls_to_fit = [mean(wls(1:10)) mean(wls(end-10:end))];
Ts1550_to_fit = [mean(Ts1550(1:10)) mean(Ts1550(end-10:end))];

bgFitP = polyfit(wls_to_fit, Ts1550_to_fit, 1);
bgFit = wls.*bgFitP(1) + bgFitP(2);

figure(); 
subplot(211);
plot(wls, Ts1550); hold on;
plot(wls, bgFit);

% Correct for bg
Ts1550 = Ts1550 ./ bgFit;
% Ts1550 = Ts1550/max(Ts1550);
subplot(212);
plot(wls, Ts1550);
xline(center_wl);

if kk >= 28 && kk <=38
    lines_detuning   = allFitParameters(kk).NewParameters(7);
%     wavelength_a     = allFitParameters(kk).NewParameters(1); % center_wl; %lambda_c; % 1.549414e-06; % lambda_c;
%     Qa               = allFitParameters(kk).NewParameters(2); %0.9e6; % Q_a;
%     Qb               = allFitParameters(kk).NewParameters(3);
%     g                = allFitParameters(kk).NewParameters(4);
%     coupling_ratio_a = allFitParameters(kk).NewParameters(5);       
%     coupling_ratio_b = allFitParameters(kk).NewParameters(6); % doesn't impact transmission lieshape
    
    wavelength_a     = center_wl; % center_wl; %lambda_c; % 1.549414e-06; % lambda_c;
    Qa               = Q_a; %mean_Qa; %0.77e6; % Q_a;
    Qb               = Q_b;
    g                = g_experiment/2/pi;
    coupling_ratio_a = Q_a/Q_a_e;% 2e6;       
    coupling_ratio_b = kappa_b_e/kappa_b;       % doesn't impact transmission lieshape
else 
    lines_detuning   = 0;
    wavelength_a     = center_wl; % center_wl; %lambda_c; % 1.549414e-06; % lambda_c;
    Qa               = Q_a; %mean_Qa; %0.77e6; % Q_a;
    Qb               = Q_b;
    g                = g_experiment/2/pi;
    coupling_ratio_a = Q_a/Q_a_e;% 2e6;       
    coupling_ratio_b = kappa_b_e/kappa_b;       % doesn't impact transmission lieshape
end

FitParameters = [wavelength_a Qa Qb  ...
                g coupling_ratio_a coupling_ratio_b ...
                lines_detuning P_pump];

Parameters = FitParameters;
Parameters(1) = center_wl;
Parameters(end) = P_pump;

%
wl_spacing      = wl_window*1.5;
freq_spacing    = 0.7e10;

center_freq = 2*pi*c./center_wl;
omegas = 2*pi*c./wls;
detuning = omegas - center_freq;

figure(all_lineshapes);
subplot(211);
plot( detuning + (iteration-1)*freq_spacing, Ts1550, '.' ); hold on; hold on;

Ts = Ts/max(smooth(Ts));
subplot(212);
plot( detuning + (iteration-1)*freq_spacing, Ts, '.' ); hold on; hold on;

if kk<38
    
    wls_sparse  = wls(1:1:end);
    omegas_sparse = 2*pi*c./wls_sparse;
    detuning_sparse = omegas_sparse - center_freq;

    SHG_fitted_fun = @(x) SHG_amplitudeNorm(Parameters, x);
    SHG_fitted = SHG_fitted_fun(wls_sparse);

    Ts_fitted_fun = @(x) SHG_pump_transmission(Parameters, x);
    Ts_fitted = Ts_fitted_fun(wls_sparse);
    
    
    subplot(211);
    plot( detuning_sparse + (iteration-1)*freq_spacing, Ts_fitted, '-', 'color', 'red');
    subplot(212);
    plot( detuning_sparse + (iteration-1)*freq_spacing, SHG_fitted, '-', 'color', 'red');
    
    plot_data{iteration}.Ts1550_fitted  = Ts_fitted;
    plot_data{iteration}.SHG_fitted     = SHG_fitted;
   
else
    
    plot_data{iteration}.Ts1550_fitted  = [];
    plot_data{iteration}.SHG_fitted     = [];
   
end

plot_data{iteration}.wls            = wls;
plot_data{iteration}.SHG            = Ts;
plot_data{iteration}.Ts1550         = Ts1550;
plot_data{iteration}.detuning       = detuning;
plot_data{iteration}.freq_spacing   = freq_spacing;

iteration = iteration + 1;

end

save(['./AllData2/lineshape_plot_data.mat'], 'plot_data');
saveas(all_lineshapes, ['./AllData2/lineshapes_model.fig']);

%% Plot theory lineshapes on top of the experimental data (actual fitted lineshapes)

all_lineshapes = figure();
iteration = 1;

modeParameters;

for kk = 1:1:length(data)

% Load single data

data_loc  = data{kk};
cal_loc   = cal{kk};

%
fiber_to_chip780 = 0.09;
fiber_to_chip1550 = 0.26;

Ts  = (data_loc.Vs(:,1) - cal_loc.backgroundV_780) * cal_loc.V_to_uW_780_output_in_fiber / fiber_to_chip780 * 1e-6;
wls = data_loc.wls*1e-9;
P_pump = data_loc.input_PM_1550 * cal_loc.W_to_uW_1550_input_in_fiber * fiber_to_chip1550 * 1e-6;
Ts1550 = (data_loc.Vs(:,2) - cal_loc.backgroundV) * cal_loc.V_to_uW_1550_output_in_fiber / fiber_to_chip1550 * 1e-6; 
TsMZI  = data_loc.Vs(:,3);

wl_filter = wls > 1549.35*1e-9 & wls < 1549.47*1e-9;    % For power 35

Ts = Ts(wl_filter);
wls = wls(wl_filter);
Ts1550 = Ts1550(wl_filter)';
TsMZI = TsMZI(wl_filter);

% MZI

MZI1550_interp   = interp1(1:length(wls), TsMZI, 1:0.1:length(wls));
wls_interp      = interp1(1:length(wls), wls, 1:0.1:length(wls));

[wls_out Ts780_norm_out delta_lambda_calibration] = mzcalibration_tmp( wls_interp*1e9, MZI1550_interp, MZI1550_interp', 0.0, 'MARKO');

wls         = interp1(1:length(wls_out), wls_out, 1:10:length(wls_out))*1e-9;
Ts1550      = Ts1550/max(Ts1550);

% Find center wavelength
[sh_max, max_ind] = max(smooth(1-Ts1550));
center_wl = mean(wls((1-Ts1550)/max(1-Ts1550) > 0.7));

% Narrow filtering
wl_window = 0.02*1e-9;
if kk >= 28 && kk <=38
    wl_filter = (wls > (center_wl-wl_window/2-0.005*1e-9)) & (wls < (center_wl+wl_window/2+0.000*1e-9));
elseif kk > 38
    wl_filter = (wls > (center_wl-wl_window/2-0.01*1e-9)) & (wls < (center_wl+wl_window/2+0.01*1e-9));
else
    wl_filter = (wls > (center_wl-wl_window/2-0.00*1e-9)) & (wls < (center_wl+wl_window/2+0.000*1e-9));
end

Ts = Ts(wl_filter);
wls = wls(wl_filter);
Ts1550 = Ts1550(wl_filter);

wls_to_fit = [mean(wls(1:10)) mean(wls(end-10:end))];
Ts1550_to_fit = [mean(Ts1550(1:10)) mean(Ts1550(end-10:end))];

bgFitP = polyfit(wls_to_fit, Ts1550_to_fit, 1);
bgFit = wls.*bgFitP(1) + bgFitP(2);

figure(); 
subplot(211);
plot(wls, Ts1550); hold on;
plot(wls, bgFit);

% Correct for bg
Ts1550 = Ts1550 ./ bgFit;
% Ts1550 = Ts1550/max(Ts1550);
subplot(212);
plot(wls, Ts1550);
xline(center_wl);

fit_threshold = 28;
if kk >= fit_threshold && kk <=38
    lines_detuning   = allFitParameters(kk).NewParameters(7);
    wavelength_a     = allFitParameters(kk).NewParameters(1); % center_wl; %lambda_c; % 1.549414e-06; % lambda_c;
    Qa               = allFitParameters(kk).NewParameters(2); %0.9e6; % Q_a;
    Qb               = allFitParameters(kk).NewParameters(3);
    g                = allFitParameters(kk).NewParameters(4);
    coupling_ratio_a = allFitParameters(kk).NewParameters(5);       
    coupling_ratio_b = allFitParameters(kk).NewParameters(6); % doesn't impact transmission lieshape
    
%     wavelength_a     = center_wl; % center_wl; %lambda_c; % 1.549414e-06; % lambda_c;
%     Qa               = Q_a; %mean_Qa; %0.77e6; % Q_a;
%     Qb               = Q_b;
%     g                = g_experiment/2/pi;
%     coupling_ratio_a = Q_a/Q_a_e;% 2e6;       
%     coupling_ratio_b = kappa_b_e/kappa_b;       % doesn't impact transmission lieshape
else 
    lines_detuning   = 0;
    wavelength_a     = center_wl; % center_wl; %lambda_c; % 1.549414e-06; % lambda_c;
    Qa               = Q_a; %mean_Qa; %0.77e6; % Q_a;
    Qb               = Q_b;
    g                = g_experiment/2/pi;
    coupling_ratio_a = Q_a/Q_a_e;% 2e6;       
    coupling_ratio_b = kappa_b_e/kappa_b;       % doesn't impact transmission lieshape
end

FitParameters = [wavelength_a Qa Qb  ...
                g coupling_ratio_a coupling_ratio_b ...
                lines_detuning P_pump];

Parameters = FitParameters;
Parameters(1) = center_wl;
Parameters(end) = P_pump;

%
wl_spacing      = wl_window*1.5;
freq_spacing    = 0.7e10;

center_freq = 2*pi*c./center_wl;
omegas = 2*pi*c./wls;
detuning = omegas - center_freq;

figure(all_lineshapes);
subplot(211);
plot( detuning + (iteration-1)*freq_spacing, Ts1550, '.' ); hold on; hold on;

Ts = Ts/max(smooth(Ts));
subplot(212);
plot( detuning + (iteration-1)*freq_spacing, Ts, '.' ); hold on; hold on;

if kk<38
    
    wls_sparse  = wls(1:1:end);
    omegas_sparse = 2*pi*c./wls_sparse;
    detuning_sparse = omegas_sparse - center_freq;

    SHG_fitted_fun = @(x) SHG_amplitudeNorm(Parameters, x);
    SHG_fitted = SHG_fitted_fun(wls_sparse);

    Ts_fitted_fun = @(x) SHG_pump_transmission(Parameters, x);
    Ts_fitted = Ts_fitted_fun(wls_sparse);
    
    
    subplot(211);
    plot( detuning_sparse + (iteration-1)*freq_spacing, Ts_fitted, '-', 'color', 'red');
    subplot(212);
    plot( detuning_sparse + (iteration-1)*freq_spacing, SHG_fitted, '-', 'color', 'red');
    
    plot_data{iteration}.Ts1550_fitted  = Ts_fitted;
    plot_data{iteration}.SHG_fitted     = SHG_fitted;
   
else
    
    plot_data{iteration}.Ts1550_fitted  = [];
    plot_data{iteration}.SHG_fitted     = [];
   
end

plot_data{iteration}.wls            = wls;
plot_data{iteration}.SHG            = Ts;
plot_data{iteration}.Ts1550         = Ts1550;
plot_data{iteration}.detuning       = detuning;
plot_data{iteration}.freq_spacing   = freq_spacing;

iteration = iteration + 1;

end

save(['./AllData2/lineshape_plot_data_fits.mat'], 'plot_data', 'fit_threshold');
saveas(all_lineshapes, ['./AllData2/lineshapes_fits.fig']);
