function [F] = fitOpticalQ(wls,Ts,wl_min,wl_max, bg_mode, do_plot)
% do a windowed fit of a lorentzian with linear background using fitlorentzian
% compute results and add to fit structure output by fitlorentzian

% window data
inds = ((wls<wl_max) & (wls> wl_min));
step = wls(10) - wls(9);
index_offset = find(abs(wls - wl_min) < step/2, 1);
wls = wls(inds);
Ts = Ts(inds);

% plot data
% h_plot = figure;


% Fit lorentzian with linear backrgound.
try
    F = fitlorentzian(wls,Ts,[],bg_mode);
catch
    error('fitlorentzian:failed_fit','Failed to perform fit.');
end



% Compute results and store them in fit structure F.
c0 = 299792458; % m/s
F.lambda_o = F.pp(2); % in nm
F.omega_o = c0/(F.pp(2)*1e-9)*2*pi;
F.Q = F.pp(2)/F.pp(3);               % lambda / FWHM

%F.peak_index = find(abs((wls - F.pp(2))) < step/2, 1); %index resonance
[~, F.peak_index] = min(abs(wls-F.pp(2)));
F.abs_peak_index = F.peak_index;

%F.abs_peak_index = F.peak_index + index_offset;

% Add Q's to title.
% title(['Q = ' num2str(F.Q) ', Qi = ' num2str(F.Q_i)])


% Plot fit and background.
    if do_plot
      figure;
      plot(wls,Ts,'.','Color',[0.4 0.4 0.9]);
      hold all
      plot(F.x,F.yfit,'-k','linewidth',1.5, 'Color', 'r');
      plot(F.x,F.bgfit,'-k','linewidth',1.5);
      xlabel('Wavelength (nm)')
      ylabel('Transmission (a.u.)')
      title(['Lambda = ' num2str(F.lambda_o) ', Q = ' num2str(round(F.Q))])
    end

end
