c = 3e8;
hbar = 1.05e-34;
omega_a = 2*pi*3e8/1549.4e-9;
omega_b = 2*omega_a;

%1550 mode, 
% Q_a = 0.605e6;
% Q_a_e = 1.001e6;
% Q_a_i = 1.529e6;
% kappa_a = omega_a/Q_a;
% kappa_a_e = omega_a/Q_a_e;

% 2021-01-18 1:16 comment out
Q_a = 0.74e6;
Q_a_i = 1.2e6;
% Q_a_e = 1.84e6;
Q_a_e = 1/(1/Q_a - 1/Q_a_i);
kappa_a = omega_a/Q_a;
kappa_a_e = omega_a/Q_a_e;

% Wider mode
% 2021-01-18 1:16 comment out
% Q_b = 0.66e6;
% Q_b_i = 0.86e6;
% Q_b_e = 2.8e6;
% Higher Q mode
% Q_b = 0.83e6;
% Q_b_i = 1.2e6;
% Q_b_e = 2.97e6;
% Q_b = 0.73e6;
% Q_b_e = Q_b/0.28;
Q_b = 0.82e6;       % 0.8e6 from fitting nonlinear lineshapes
Q_b_i = 1.2e6;
Q_b_e = 1/(1/Q_b - 1/Q_b_i);

kappa_b = omega_b/Q_b;
kappa_b_e = omega_b/Q_b_e;

%1e-6;
duty_cycle = 0.76;
g_theory = 2*pi*250e3 * sin(pi*duty_cycle);
% g_theory = 2*pi*180e3*0.8;

% 2021-01-18 1:16 comment out
% g_experiment = 2*pi* 1.0905e+05; % from slope efficiency
g_experiment = 2*pi* 1.3e+05; % from slope efficiency
