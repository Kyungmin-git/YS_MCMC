function [G, F, Cdisk] = compute_G_hankel(r_obs, R, H, D, E, nu)
% compute_G_hankel
%
% Computes the vertical surface-displacement transfer function G
% for an equivalent axisymmetric circular source whose total
% pressure-induced volume change is constrained by the oblate
% spheroidal cavity solution.
%
% u_z(omega) = G * DeltaP(omega)
%
% INPUTS
%   r_obs : horizontal radial offset from source center [m]
%   R     : horizontal semi-axis / equivalent disk radius [m]
%   H     : source center depth [m]
%   D     : full gas-pocket thickness [m]
%           (oblate vertical semi-axis c = D/2)
%   E     : Young's modulus [Pa]
%   nu    : Poisson's ratio
%
% OUTPUTS
%   G     : pressure-to-vertical-displacement transfer function [m/Pa]
%   F     : dimensionless geometric Hankel response
%   Cdisk : effective distributed-disk specific compliance [Pa^-1]

    %% ----------------------------------------------------
    % Elastic constants
    % -----------------------------------------------------

    mu = E / (2 * (1 + nu));

    %% ----------------------------------------------------
    % Oblate-derived equivalent disk compliance
    %
    % Oblate spheroid:
    %   a = b = R
    %   c = D/2
    %
    % Cdisk = (2/3) Cobl
    % so that the equivalent uniform disk preserves the
    % total pressure-induced volume change DeltaV.
    % -----------------------------------------------------

    Cdisk = (1 / mu) * ...
            ( (8 * (1 - nu) * R) / (3 * pi * D) ...
              - (1 - 2*nu) / (1 + nu) );

    %% ----------------------------------------------------
    % Dimensionless geometry
    % -----------------------------------------------------

    s = r_obs / R;
    h = H / R;

    %% ----------------------------------------------------
    % Hankel geometric response
    %
    % F(s,h) = integral_0^Inf
    %          exp(-h*q) J1(q) J0(s*q) dq
    % -----------------------------------------------------

    if abs(s) < 1e-12

        % Exact centered solution
        %
        % F(0,h) = 1 - h/sqrt(1+h^2)

        F = 1 - h / sqrt(1 + h^2);

    else

        integrand = @(q) exp(-h .* q) .* ...
                         besselj(1, q) .* ...
                         besselj(0, s .* q);

        F = integral(integrand, 0, Inf, ...
                     'RelTol', 1e-7, ...
                     'AbsTol', 1e-10);

    end

    %% ----------------------------------------------------
    % Pressure-to-displacement transfer function
    %
    % u_z(omega) = G * DeltaP(omega)
    %
    % G = 2(1-nu) Cdisk D F
    % -----------------------------------------------------

    G = 2 * (1 - nu) * Cdisk * D * F;

end