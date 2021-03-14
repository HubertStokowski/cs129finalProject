function [Wave_out Signal_out delta_lambda_calibration] = mzcalibration_tmp(Wave_in,Signal_in,Calib_in,offset,which_MZI)
% Qiang     10/13/2010, code modified from Chris's code which has an incorrect background removal for MZI interferences.
% Function for calibrating the wavelength data using the MZ interferometer.
% Usage: [Wave_out Signal_out] = MZIWaveCalib(Wave_in,Signal_in,Calib_in,offset)
%           Wave_in   = wavelength data vector  (in NANOMETERS!)
%           Signal_in = transmission/reflection data
%           Calib_in  = MZ interference fringes data
%           offset    = percentage of data at beginning and end (0-50%) of 
%                       trace to be dropped (motor acceleration, etc.)
% ... MZI FSR is calibrated by network analyzer at both 980nm and 1530nm bands. 

c=3*10^11;              % velocity of light (nm/us)

output_graph = 0; % 1 for output graph

%% Calibrated free-spectral ranges. The two numbers are quite close because of the flat group index of silica fiber over broad spectral range. 
if strcmpi(which_MZI,'BLUE')
FSR_980=220.22;                                             % (MHz) FSR at 979.5nm, calibrated by network analyzer Agilent 4396B.
FSR_1530=220.31;                                            % (MHz) FSR at 1531.0nm, calibrated by network analyzer Aiglent 4396B. 
FSR_1300 = mean([FSR_1530 FSR_980]);%221.96;                                          % (MHZ) FSR at 1320nm, calibrated using wavemeter + fringe counting
elseif strcmpi(which_MZI,'GHENT')
FSR_1530 = 453.6;%(MHz) calibrated by Tim Blasius 
FSR_980 = 453.6; %(MHz) guessed by Alex Krause based on previous 1530->980 conversion values. 
FSR_1300 = mean([FSR_1530 FSR_980]);
elseif strcmpi(which_MZI,'ASN_MZI')
FSR_1530 = 459.5; %453.6; %(MHz) calibrated by Tim Blasius 
FSR_980 = 459.5; %453.6; %(MHz) guessed by Alex Krause based on previous 1530->980 conversion values. 
FSR_1300 = mean([FSR_1530 FSR_980]);
elseif strcmpi(which_MZI,'BIGWHITE')
FSR_1530 = 325.28; %(MHz) calibrated by Tim Blasius 
FSR_980 = 325.18; %(MHz) guessed by Alex Krause based on previous 1530->980 conversion values. 
FSR_1300 = mean([FSR_1530 FSR_980]);
elseif strcmpi(which_MZI,'SMALLWHITE')
FSR_1530 = 173; %(MHz) 
FSR_980 = 173; %(MHz)
FSR_1300 = mean([FSR_1530 FSR_980]);
elseif strcmpi(which_MZI,'SNOWWHITE')
FSR_1530 = 1.130784270099993e+002; %(MHz) 
FSR_980 = 1.130784270099993e+002; %(MHz) 
FSR_1300 = mean([FSR_1530 FSR_980]);
elseif strcmpi(which_MZI,'MARKO')
FSR_1530 = 67.69103657394908; %(MHz) 
FSR_980 = 67.69103657394908; %(MHz) 
FSR_1300 = mean([FSR_1530 FSR_980]);
elseif strcmpi(which_MZI,'780')
FSR_1530 = 39.9099; %(MHz) 
FSR_980 = 39.9099; %(MHz) 
FSR_1300 = mean([FSR_1530 FSR_980]);
end
% FSR=0.001*1.57; % (1.57 +/- 0.03) pm @ 1460.5nm is the calibration on Matt's L3 pg.143

OrderPoly_BackgroundFit=6;                                  % polynomial order of the background removal for MZI interference

%% Chop data if desired
if offset==0
    Wave1 = Wave_in;
    Signal1 = Signal_in;
    Calib1 = Calib_in;
else
    cut = floor(length(Wave_in)*(offset/100));
    Wave1 = Wave_in(cut:(end-cut));
    Signal1 = Signal_in(cut:(end-cut));
    Calib1 = Calib_in(cut:(end-cut));
end

WavelengthIndex=[1:length(Calib1)]';
MidIndex0=round(length(WavelengthIndex)/2);
WaveIndex0=WavelengthIndex(MidIndex0);
Wavelength0=Wave1(MidIndex0);                              % center wavelength of the piezo scanning range (nm), which is used as the reference wavelength.
Frequency0=c/Wavelength0;                                  % center frequency (MHz)

if Wavelength0 < 1250                                   % if the center wavelength is smaller than 1300nm, use the calibrated FSR at 980nm band  
    FSR=FSR_980;
elseif Wavelength0 >= 1350                              % if the center wavelength is larger than 1300nm, use the calibrated FSR at 1530nm band
    FSR=FSR_1530;
else
    FSR=FSR_1300;
end

%% Remove the background on the MZI interference
CropPercentage=0.04;        % percentage of the total piezo scanning range at the scanning edge to be excluded from finding the MZI background, to prevent errors in finding fringe peaks/valleys.
IndexCrop=find(WavelengthIndex >= CropPercentage*length(Calib1) & WavelengthIndex <= (1-CropPercentage)*length(Calib1));        % find the interference peaks and valleys with 5-95% scanning range
WaveIndexCrop=WavelengthIndex(IndexCrop);
CalibCrop=Calib1(IndexCrop);
CalibCropSmooth=smooth(smooth(CalibCrop));                  % smooth the interference pattern to remove some noises
[MZIFringePeaks PeakLocs]=findpeaks(CalibCropSmooth);       % search for the transmission peaks of MZI interferences
WaveIndex_MZIPeak=WaveIndexCrop(PeakLocs);
[MMID ValleyLocs]=findpeaks(-CalibCropSmooth);              % search for the transmission valley of MZI interferences
WaveIndex_MZIValley=WaveIndexCrop(ValleyLocs);
MZIFringeValleys=CalibCropSmooth(ValleyLocs);

% find the transmission background along the interference peaks, fitted with polynomial
MidIndex_MZIPeak=round(length(WaveIndex_MZIPeak)/2);
WaveIndex_MZIPeak0=WaveIndex_MZIPeak(MidIndex_MZIPeak);
DWaveIndex_MZIPeak=WaveIndex_MZIPeak-WaveIndex_MZIPeak0;
[BackgroundPeak_fit,STE_Peak,muTE_Peak]=polyfit(DWaveIndex_MZIPeak(:),MZIFringePeaks(:),OrderPoly_BackgroundFit);
MZIBackgroundPeak_Fit=polyval(BackgroundPeak_fit,WavelengthIndex-WaveIndex_MZIPeak0,[],muTE_Peak);

% find the transmission background along the interference valleys, fitted with polynomial
MidIndex_MZIValley=round(length(WaveIndex_MZIValley)/2);
WaveIndex_MZIValley0=WaveIndex_MZIValley(MidIndex_MZIValley);
DWaveIndex_MZIValley=WaveIndex_MZIValley-WaveIndex_MZIValley0;
[BackgroundValley_fit,STE_Valley,muTE_Valley]=polyfit(DWaveIndex_MZIValley(:),MZIFringeValleys(:),OrderPoly_BackgroundFit(:));
MZIBackgroundValley_Fit=polyval(BackgroundValley_fit,WavelengthIndex-WaveIndex_MZIValley0,[],muTE_Valley);

MZIBackground=(MZIBackgroundPeak_Fit+MZIBackgroundValley_Fit)/2;        % MZI background defined as the mean between the fitted peaks and valleys
MZIFringeAmplitude=(MZIBackgroundPeak_Fit-MZIBackgroundValley_Fit)/2;   % magnitudes of MZI interference fringes defined as half of the difference between the fitted peaks and valleys

% normalize MZI transmission and bring down to 0 mean.
% Calib_norm=Calib1./MZIBackground-1;                                   % normalize by directly divided by the mean background. 
Calib_norm=(Calib1-MZIBackground)./MZIFringeAmplitude;                  % normalize by subtracted by the mean background, and divided by the fringe amplitudes.

Calib_norm=smooth(Calib_norm);                                          % smooth the interference pattern to remove some noises

% figure(1);
% plot(WavelengthIndex,calib1,'b',WaveIndex_MZIPeak,MZIFringePeaks,'r',WaveIndex_MZIValley,MZIFringeValleys,'c',WavelengthIndex,MZIBackground,'k');
% figure(2);
% plot(WavelengthIndex,Calib_norm);


%% Find zeros and interpolate wavelengths
% find the zeros of the fringes, interpolated as the points where the crossing occurs
Zeros_ind = [];
for m=2:length(Wave1)
    if Calib_norm(m-1)*Calib_norm(m) < 0
        Zeros_ind(end+1,:) = [interp1(Calib_norm((m-1):m),[(m-1) m],0) 0];
    end
end

LocationVector=[1:size(Zeros_ind,1)]';
CenterWavelengthLocation=interp1(Zeros_ind(:,1),LocationVector,WaveIndex0);     % find the relative locations of the reference center wavelength in the calibrated fringe zeros.

Zeros_ind(:,2)=Frequency0-0.5*FSR*(LocationVector-CenterWavelengthLocation);    % calibrated frequencis of MZI interference zeros, used center frequency as the reference frequency.
Zeros_ind(:,2)=c./Zeros_ind(:,2);                                               % change unit to be wavelnegth (nm)

Wave_out = interp1(Zeros_ind(:,1),Zeros_ind(:,2),1:length(Wave1),'spline','extrap')';   % calibrated wavelength (nm)
Signal_out = Signal1;                                                                   % the signal remains intact.

delta_lambda_calibration = Wave_out-Wave1;
%% output graph
if output_graph ==1
%     SetFigPos(2,1234)
%     clf(1);
    subplot(4,1,1:3)
    box on; % hold on;
    plot(Wave_in,Signal_in,'b.-',Wave_out,Signal_out,'m.-',Wave_out,Calib_norm,'g.-',Zeros_ind(:,2),0,'go');
    legend('original','calibrated','fringes');
    legend('Location','SouthEast')
    xlabel('Wavelength (nm)')
    line([(min([Wave_in(1) Wave_out(1)])-10) (max([Wave_in(end) Wave_out(end)])+10)],...
            [0 0],'Color',[0 0 0]);
%     plot(Zeros_ind(:,2),0,'go');
    xlim([min([Wave_in(1) Wave_out(1)]) max([Wave_in(end) Wave_out(end)])])
    %hold off
    subplot(4,1,4)
    plot(Wave_out,Wave_out-Wave1,'b.-');
    xlabel('Calibrated Wavelength (nm)');
    ylabel('\lambda_c_a_l_i_b - \lambda_d_a_t_a (nm)')
    grid on;
    xlim([min([Wave_in(1) Wave_out(1)]) max([Wave_in(end) Wave_out(end)])])
end