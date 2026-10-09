function [Posc, uspect, uspectf, anglef, synthetic, A_res, A_exc, ...
    A_path, synangle, G] = forwardmodelgaspocket(m, nu, YE, t0, qn, sigma_k_log)
%FORWARDMODELGASPOCKET Yellowstone gas-pocket pressure/displacement model.

%   Inputs:
%     m = [T, Q, R, L, -log10(kappa), -log10(D), log10(phi_percent), r_obs]
%     T: K; Q: kg/s; R, L, D, r_obs: m; kappa: m^2.
%     nu: Poisson ratio; YE: Young's modulus (Pa).
%     t0: excitation times (s); qn: corresponding excitation masses (kg).
%     sigma_k_log: log10 of the excitation-width coefficient.
%
%   Outputs retain the original order and vector orientations:
%     Posc: pressure relative to equilibrium (Pa), row vector.
%     uspect: unfiltered displacement-spectrum magnitude, row vector.
%     uspectf: filtered positive-frequency magnitude, column vector.
%     anglef: full filtered FFT phase (rad), column vector.
%     synthetic: unfiltered displacement waveform (m), row vector.
%     A_res, A_exc: response and excitation spectra, row vectors.
%     A_path: zero-filled compatibility output (not used in this model).
%     synangle: unfiltered displacement-spectrum phase, row vector.
%     G: elastic pressure-to-displacement transfer coefficient.


%% Physical parameters
mu_g = 1.5e-5; % Gas viscosity (Pa s).
M = 0.018;     % Water-vapor molar mass (kg/mol).
Rg = 8.3145;   % Universal gas constant (J/(mol K)).
Pex = 1.5e6;   % External pressure (Pa).

T = m(1);
Q = m(2);
R = m(3);
L = m(4);
kappa = 10^(-m(5));
D = 10^(-m(6));
phi = 0.01*10^(m(7));
r_obs = m(8);

S = pi*R^2;

%% Auxiliary parameters
beta_a = S*phi*M/(Rg*T);
beta_b = mu_g*phi/(kappa*(Pex-Rg*T*Q^2/(S^2*phi^2*M*Pex)));
beta_c = Pex*M/(Rg*T*(Pex-Rg*T*Q^2/(S^2*phi^2*M*Pex)));
beta_d = 2*Q/(S*phi*(Pex-Rg*T*Q^2/(S^2*phi^2*M*Pex)));
beta_e = S*M*D/(Rg*T);
P0 = Pex+mu_g*Rg*T*Q*L/(S*kappa*M*(Pex-Rg*T*Q^2/(S^2*phi^2*M*Pex)));

%% Harmonic-oscillator coefficients
GAMMA0 = 1;
GAMMA1 = (2*(beta_a*beta_d+beta_b*beta_e)*L+beta_a*beta_b*L^2)/(2*beta_a);
GAMMA2 = (2*beta_c*beta_e*L+beta_a*beta_c*L^2)/(2*beta_a);
gamma0 = beta_b*L/beta_a;
gamma1 = beta_c*L/beta_a;

%% Simulation grid (80 Hz)
global LPstartt LPendt tau

tau = LPendt-LPstartt;
sampling_frequency = 80;
dt = 1/sampling_frequency;
time1 = 0:dt:tau;
Lw = length(time1);
df = (1/dt)/Lw;

%% Elastic surface-displacement transfer
% Axisymmetric circular-source Hankel/Bessel formulation: u_z = G*DeltaP.
[G, ~, ~] = compute_G_hankel(r_obs, R, L, D, YE, nu);

%% Pressure and displacement spectra
u_z = zeros(1, Lw);
A_p = zeros(1, Lw);
A_res = zeros(1, Lw);
A_exc = zeros(1, Lw);
A_path = zeros(1, Lw); % Retained for output-interface compatibility.

sigma_k = 10^(sigma_k_log);
sigma_row = sigma_k*qn.^(1/3);

% Retain the original 1/tau spacing and zero-filled negative-frequency
% entries for consistency with the inversion's existing numerical model.
mm = 1;
for frequency = 0:1/tau:1/(2*dt)
    omega = 2*pi*frequency;

    A_res(mm) = ...
        (gamma0*(1i*omega)^0+gamma1*(1i*omega)^1) / ...
        (GAMMA0*(1i*omega)^0+GAMMA1*(1i*omega)^1+GAMMA2*(1i*omega)^2);

    A_exc(mm) = sum(qn .* exp(-1i*omega*t0) .* ...
        exp(-(1/2)*sigma_row.^2.*omega^2));

    A_p(mm) = A_res(mm)*A_exc(mm);
    u_z(mm) = G*A_p(mm);
    mm = mm+1;
end

%% Time-domain reconstruction (original normalization)
P = ifft(A_p, 'symmetric')*sqrt(Lw)/sqrt(dt/df);
Posc = P+Pex-P0;
uspect = abs(u_z);
U_z = ifft(u_z, 'symmetric')*sqrt(Lw)/sqrt(dt/df);
synangle = angle(u_z);
synthetic = real(U_z);

%% Filtered displacement spectrum (1–35 Hz)
order = 4;
low_frequency = 1;
high_frequency = 35;
[b, a] = butter(order/2, ...
    [low_frequency high_frequency]/(sampling_frequency/2));
syntheticf = filtfilt(b, a, synthetic);
Lw = length(syntheticf);
NFFT = Lw;
water = fft(syntheticf', NFFT)/sqrt(Lw)*sqrt(dt/df);
% floor explicitly expresses the original colon's last positive-bin index.
uspectf = abs(water(1:floor(NFFT/2)+1));
anglef = angle(water);
end
