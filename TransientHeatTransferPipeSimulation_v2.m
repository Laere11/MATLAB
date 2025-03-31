% TransientHeatTransferPipeSim.m
% Simulates 1D transient heat conduction in a pipe with partial cooling.
% Modified parameters to enhance transient dynamics.
clear; clc; close all;

%% USER INPUTS
total_time = 500;        % Total simulation time in seconds
frame_interval = 5;      % Plot update interval in seconds

%% Physical and Material Properties
L = 3;                  % Length of the pipe (m)
D = 0.05;                % Diameter of the pipe (m)

% Using aluminum-like properties for faster conduction:
k = 1205;                 % Thermal conductivity (W/m·K)
rho = 2700;              % Density (kg/m³)
cp = 900;                % Specific heat (J/kg·K)

% Adjusted cooling conditions for a more gradual change:
T_initial = 40;          % Initial temperature (°C)
T_c = 25;                % Cooling region temperature (°C)
fractionCool = 0.01;      % Only 10% of the pipe (near inlet) is cooled

%% Numerical Parameters
N = 101;                         % Number of spatial nodes
dx = L / (N - 1);                % Spatial resolution
x = linspace(0, L, N)';          % Spatial grid

alpha = k / (rho * cp);          % Thermal diffusivity

% Stability condition for explicit FTCS
dt_stable = dx^2 / (2 * alpha);
safety_factor = 0.4;
dt_calc = safety_factor * dt_stable;

% Manually override dt for better animation control:
manual_dt = 1; % Use 1 second per time step
dt = min(dt_calc, manual_dt);

n_steps = max(round(total_time / dt), 1);
frame_step = max(round(frame_interval / dt), 1);

%% Initialization
T = T_initial * ones(N,1);       % Set initial temperature across the pipe
cool_nodes = round(fractionCool * N);
T(1:cool_nodes) = T_c;           % Apply cooling to the inlet region

%% Time-Stepping Loop (Explicit FTCS)
figure;
for step = 1:n_steps
    T_old = T;
    
    % Update interior nodes using FTCS scheme:
    for i = 2:N-1
        T(i) = T_old(i) + alpha * dt/(dx^2) * (T_old(i+1) - 2*T_old(i) + T_old(i-1));
    end
    
    % Re-apply boundary conditions:
    T(1:cool_nodes) = T_c;    % Cooling at inlet remains fixed
    T(end) = T_old(end);      % Insulated outlet (Neumann BC)
    
    % Plot updates at specified intervals:
    if mod(step, frame_step) == 0
        plot(x, T, 'LineWidth', 2);
        xlabel('Pipe Length (m)');
        ylabel('Temperature (°C)');
        title(sprintf('Time = %.1f sec', step * dt));
        grid on;
        drawnow;
    end
end
