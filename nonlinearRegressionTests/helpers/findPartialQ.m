function [Qtot, Qi, Qe] = findPartialQ(F, deviceType, varargin)
% %3/29/2018 TPM modified to work with uwave mode object from
% %fitMicrowaveQ()

%10/11/2017 TPM
%F is the struct returned by fitopticalQ()
%deviceType can be 'symmRing' or 'symmCavityTrans' or 'symmCavityRefl' or 'singleSidedCavity'
%MUST supply grating_profile for 'symmCavityTrans'

if (strcmp(deviceType,'symmCavityTrans') & length(varargin) ~= 1)
  error('You must supply the grating_profile for device type symmCavityTrans');
end

Qtot = F.Q;
omega_o = F.omega_o;
peak_index = F.peak_index;
abs_peak_index = F.abs_peak_index;
%c = 299792458; % m/s
k_tot = omega_o/Qtot; %kappa total in radians/s

switch deviceType
  case 'symmCavityTrans'
    %symmetric cavity, transmission, any coupling
    grating_profile = varargin{1};
    T_o = max(F.yfit) / grating_profile(abs_peak_index);
    k_e = sqrt(T_o) * k_tot / 2;
    k_i = k_tot - 2*k_e;
  case 'symmCavityRefl'
    %symmetric cavity, reflection
    R_o = (F.bgfit(peak_index) + F.pp(1)) / F.bgfit(peak_index);
    
    %undercoupled
    k_i = k_tot * sqrt(R_o);
    k_e = (k_tot - k_i) / 2;
    
    %overcoupled
    %temp = k_i;
    %k_i = k_e;
    %k_e = temp;
        
    
  case 'singleSidedCavity'
    %infinite end mirror, reflection
    R_o = (F.bgfit(peak_index) + F.pp(1)) / F.bgfit(peak_index);
    
    %undercoupled
    k_e = (1 - sqrt(R_o)) * k_tot / 2;
    k_i = k_tot - k_e;
    
    %overcoupled
    %k_i = (1 - sqrt(R_o)) * k_tot / 2; 
    %k_e = k_tot - k_i;
  case 'symmRing'
    T_o = (F.bgfit(peak_index) + F.pp(1)) / F.bgfit(peak_index);
    k_i = k_tot * sqrt(T_o);
    k_e = (k_tot - k_i) / 2;
    
  case 'ssuc'
    T_o = (F.bgfit(peak_index) + F.pp(1)) / F.bgfit(peak_index);
    %k_i = k_tot * sqrt(T_o);
    %k_e = (k_tot - k_i) ; %don't divide by two here
    
    %over-coupled
    %k_i = k_tot*(1-sqrt(T_o))/2;
    %k_e = k_tot*(1+sqrt(T_o))/2;
    
    %under-coupled
    k_e = k_tot*(1-sqrt(T_o))/2;
    k_i = k_tot*(1+sqrt(T_o))/2;

  case 'ssoc'
    T_o = (F.bgfit(peak_index) + F.pp(1)) / F.bgfit(peak_index);
    %k_i = k_tot * sqrt(T_o);
    %k_e = (k_tot - k_i) ; %don't divide by two here
    
    %over-coupled
    k_i = k_tot*(1-sqrt(T_o))/2;
    k_e = k_tot*(1+sqrt(T_o))/2;
    
    %under-coupled
%     k_e = k_tot*(1-sqrt(T_o))/2;
%     k_i = k_tot*(1+sqrt(T_o))/2;
end

Qe = omega_o/k_e;
Qi = omega_o/k_i;

end
