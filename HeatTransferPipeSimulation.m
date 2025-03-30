% HeatTransferPipeSimulation.m
% This script simulates the temperature distribution along a pipe containing 
% liquid methane fuel, with a portion of the pipe subjected to cryogenic cooling.
% The cooling is modeled as a convective heat loss term applied in the cooled region.
% The simulation uses a finite difference method to solve the steady-state heat conduction
% equation. A parameter study is included to see how pipe diameter and length affect
% the temperature gradient.

clear; clc; close all;

%% Parameters

% Pipe geometry and properties
L = 10;            % Total pipe length in meters
D = 0.05;          % Pipe diameter in meters

% Material properties and heat transfer coefficients
k = 0.1;           % Thermal conductivity (W/m*K) (example value, adjust as needed)
h = 50;            % Convective heat transfer coefficient in the cooled region (W/m^2*K)
% Biot number for the cooled region (affects the rate of convective cooling)
Bi = 4*h/(k*D);

% Temperature conditions
T_initial = 20;    % Temperature at the pipe inlet (°C)
T_c = -80;         % Cryogenic cooling temperature (°C) applied in the cooled region

% Fraction of the pipe subjected to cryogenic cooling (starting from the inlet)
fractionCool = 0.3;  % For example, 30% of the pipe is cooled

%% Discretization of the domain
N = 101;                     % Number of discretization nodes
dx = L/(N-1);                % Spatial step size
x = linspace(0, L, N)';       % Spatial grid

%% Assemble the coefficient matrix A and the right-hand side vector b
% The steady-state heat conduction equation is discretized as:
%   For cooled region: (T(i-1) - 2*T(i) + T(i+1))/dx^2 - Bi*(T(i) - T_c) = 0
%   For uncooled region: (T(i-1) - 2*T(i) + T(i+1))/dx^2 = 0

A = zeros(N, N);
b = zeros(N, 1);

% Boundary condition at the inlet (x = 0): fixed temperature
A(1,1) = 1;
b(1) = T_initial;

% Internal nodes
for i = 2:N-1
    if x(i) <= fractionCool * L
        % Node is in the cooled region
        A(i, i-1) = 1/dx^2;
        A(i, i)   = -2/dx^2 - Bi;
        A(i, i+1) = 1/dx^2;
        b(i) = -Bi * T_c;  % Forcing term due to cooling
    else
        % Node is in the uncooled region (pure conduction)
        A(i, i-1) = 1/dx^2;
        A(i, i)   = -2/dx^2;
        A(i, i+1) = 1/dx^2;
        b(i) = 0;
    end
end

% Neumann boundary condition at the outlet (x = L): insulated end (dT/dx = 0)
% Approximated as: (T(N) - T(N-1))/dx = 0  ->  T(N) = T(N-1)
A(N, N-1) = -1/dx;
A(N, N)   = 1/dx;
b(N) = 0;

%% Solve the linear system A*T = b for the temperature distribution
T = A\b;

%% Plot the temperature distribution for the baseline case
figure;
plot(x, T, 'LineWidth', 2);
xlabel('Position along the pipe (m)');
ylabel('Temperature (°C)');
title('Temperature Distribution along the Pipe');
grid on;

%% Parameter Study: Varying Pipe Diameter and Length
% Here we simulate two different diameters and two different lengths to see
% how these parameters affect the temperature profile.

diameters = [0.05, 0.1]; % Pipe diameters in meters
lengths = [10, 20];      % Pipe lengths in meters

figure;
for i = 1:length(diameters)
    for j = 1:length(lengths)
        D_local = diameters(i);
        L_local = lengths(j);
        % Recalculate Biot number for the current diameter
        Bi_local = 4*h/(k*D_local);
        % Discretize the domain for the current length
        N_local = 101;
        dx_local = L_local/(N_local-1);
        x_local = linspace(0, L_local, N_local)';
        
        % Assemble the linear system for the current case
        A_local = zeros(N_local, N_local);
        b_local = zeros(N_local, 1);
        
        % Inlet boundary condition
        A_local(1,1) = 1;
        b_local(1) = T_initial;
        
        for ii = 2:N_local-1
            if x_local(ii) <= fractionCool * L_local
                A_local(ii, ii-1) = 1/dx_local^2;
                A_local(ii, ii)   = -2/dx_local^2 - Bi_local;
                A_local(ii, ii+1) = 1/dx_local^2;
                b_local(ii) = -Bi_local * T_c;
            else
                A_local(ii, ii-1) = 1/dx_local^2;
                A_local(ii, ii)   = -2/dx_local^2;
                A_local(ii, ii+1) = 1/dx_local^2;
                b_local(ii) = 0;
            end
        end
        
        % Outlet boundary condition: insulated (Neumann)
        A_local(N_local, N_local-1) = -1/dx_local;
        A_local(N_local, N_local)   = 1/dx_local;
        b_local(N_local) = 0;
        
        % Solve for temperature distribution
        T_local = A_local\b_local;
        
        % Plot the results in a subplot
        subplot(length(diameters), length(lengths), (i-1)*length(lengths) + j);
        plot(x_local, T_local, 'LineWidth', 2);
        xlabel('x (m)');
        ylabel('T (°C)');
        title(sprintf('D = %.2f m, L = %.1f m', D_local, L_local));
        grid on;
    end
end

% End of HeatTransferPipeSimulation.m
