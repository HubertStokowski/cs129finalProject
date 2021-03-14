function [P_SHG_depl, P_fund_t, beta_depl, alpha_depl] = SHG_amplitudeNorm(Parameters, pump_wl)
% Parameters are assigned as follows:
% Parameters(1) = wavelength_a
% Parameters(2) = Qa
% Parameters(3) = Qb
% Parameters(4) = g
% Parameters(5) = coupling_ratio_a
% Parameters(6) = coupling_ratio_b
% Parameters(7) = lines_detuning
% Parameters(8) = P_pump

    wavelength_a        = Parameters(1);
    Qa                  = Parameters(2);
    Qb                  = Parameters(3);
    g                   = Parameters(4);
    coupling_ratio_a    = Parameters(5);
    coupling_ratio_b    = Parameters(6);
    lines_detuning      = Parameters(7);
    P_pump              = Parameters(8);
    
    c               = 3e8;
    hbar            = 1.0545718e-34;

    omega_a = 2*pi*c/ wavelength_a;
    omega_b = omega_a * 2 + lines_detuning;
    omega_p = 2*pi*c./ pump_wl;
    
    g       = 2 * pi * g;

    % 1550 mode
    Q_a = Qa;
    k_a = omega_a / Q_a;
    k_a_e = coupling_ratio_a * k_a;

    % 780 mode
    Q_b = Qb;
    k_b = omega_b / Q_b;
    k_b_e = coupling_ratio_b * k_b;

    % omega_p = omega_a + linspace(-5,0,50)*k_a;

    % Delta_a = omega_a - omega_p;
    % Delta_b = omega_b - 2*omega_p;

    Delta_a = omega_a - omega_p;
    Delta_b = omega_b - 2 * omega_p;

    %% Calculate SHG power in resonator

    P_pump = P_pump;  % Watts
   
    parfor jj = 1:length(Delta_a)
        %pump field inside the ring, coupling in from P_pump
        E_pump = sqrt(P_pump * k_a_e / hbar / omega_a);
        polynomial_coeff = [ (2*g^2/(1i*Delta_b(jj) + k_b/2)) ...
                                0        ...
                                (1i*Delta_a(jj) + k_a/2)   ...
                                 1i*E_pump ];
        xx = sym('xx')
        sol(jj) = solve(polynomial_coeff(1) * abs(xx)^2*xx + polynomial_coeff(3) * xx + polynomial_coeff(4) == 0);
        alpha_depl(jj) = double(sol(jj));
        beta_depl(jj) = 1i*g*alpha_depl(jj)^2 / (-1i*Delta_b(jj) - k_b/2);
        ampl_SHG_depl(jj) = sqrt(k_b_e * hbar*omega_b) * beta_depl(jj);
        P_SHG_depl(jj) = k_b_e * abs(beta_depl(jj))^2 * hbar*omega_b;
        P_fund_t(jj)   = abs(sqrt(P_pump/hbar/omega_a) - 1i*sqrt(k_a_e) * alpha_depl(jj) ).^2 * hbar*omega_a;
        
    end

    % Normalize peak to 1
    P_SHG_depl = P_SHG_depl/max(P_SHG_depl);
    
end