% HeatTransferPipeSimulation_v2.m
% This updated version includes more realistic values for material conductivity,
% cooling strength, and the extent of the cooled region. It computes the steady-state
% temperature distribution along a cryogenically cooled pipe using finite difference methods.

clear; clc; close all;

%% Updated Parameters

% Pipe geometry and properties
L = 10;            % Total pipe length in meters
D = 0.05;          % Pipe diameter in meters

% Material properties and heat transfer coefficients (UPDATED)
k = 1.0;           % Thermal conductivity (W/m*K)
h = 10;            % Convective heat transfer coefficient (W/m^2*K)

Bi = 4*h/(k*D);    % Updated Biot number based on new k and h

% Temperature conditions
T_initial = 20;    % Temperature at pipe inlet (°C)
T_c = -80;         % Cryogenic cooling temperature (°C)

% Fraction of pipe under cryogenic cooling (UPDATED)
fractionCool = 0.2;  % 20% of the pipe is cooled

%% Discretization of the domain
N = 101;                     % Number of discretization nodes
dx = L/(N-1);                % Spatial step size
x = linspace(0, L, N)';      % Spatial grid

%% Assemble the coefficient matrix A and the right-hand side vector b
A = zeros(N, N);
b = zeros(N, 1);

% Boundary condition at the inlet (x = 0): fixed temperature
A(1,1) = 1;
b(1) = T_initial;

% Internal nodes
for i = 2:N-1
    if x(i) <= fractionCool * L
        % Cooled region
        A(i, i-1) = 1/dx^2;
        A(i, i)   = -2/dx^2 - Bi;
        A(i, i+1) = 1/dx^2;
        b(i) = -Bi * T_c;
    else
        % Uncooled region
        A(i, i-1) = 1/dx^2;
        A(i, i)   = -2/dx^2;
        A(i, i+1) = 1/dx^2;
        b(i) = 0;
    end
end

% Neumann boundary condition at the outlet (x = L): insulated end
A(N, N-1) = -1/dx;
A(N, N)   = 1/dx;
b(N) = 0;

%% Solve the linear system A*T = b
T = A\b;

%% Plot the temperature distribution for the baseline case
figure;
plot(x, T, 'b-', 'LineWidth', 2);
xlabel('Position along the pipe (m)');
ylabel('Temperature (\circC)');
title('Temperature Distribution along the Pipe (Updated)');
grid on;

%% Parameter Study: Varying Pipe Diameter and Length
diameters = [0.05, 0.1]; % Pipe diameters in meters
lengths = [10, 20];      % Pipe lengths in meters

figure;
for i = 1:length(diameters)
    for j = 1:length(lengths)
        D_local = diameters(i);
        L_local = lengths(j);
        Bi_local = 4*h/(k*D_local);
        N_local = 101;
        dx_local = L_local/(N_local-1);
        x_local = linspace(0, L_local, N_local)';

        A_local = zeros(N_local, N_local);
        b_local = zeros(N_local, 1);

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

        A_local(N_local, N_local-1) = -1/dx_local;
        A_local(N_local, N_local)   = 1/dx_local;
        b_local(N_local) = 0;

        T_local = A_local\b_local;

        subplot(length(diameters), length(lengths), (i-1)*length(lengths) + j);
        plot(x_local, T_local, 'LineWidth', 2);
        xlabel('x (m)'); ylabel('T (\circC)');
        title(sprintf('D = %.2f m, L = %.1f m', D_local, L_local));
        grid on;
    end
end
