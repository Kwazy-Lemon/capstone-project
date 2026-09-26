%% Transient Analysis using Modified Nodal Analysis (MNA)
% Week 1 - Linear Circuit Prototype
%
% Circuit:
%
%                 R1 = 1 kOhm
% Vin(t) --------/\/\/\--------o Vout
%                               |
%                               +---- R2 = 2 kOhm ---- GND
%                               |
%                               +---- C = 100 nF ----- GND
%
% Input:
%   0 V before t = 0
%   10 V step at t = 0
%
% Transient analysis is performed using the
% Backward Euler method.

clear;
clc;
close all;

%% Circuit Parameters

Vs = 10;                 % Step voltage [V]

R1 = 1e3;                % R1 [Ohm]
R2 = 2e3;                % R2 [Ohm]
C  = 100e-9;             % C [F]

G1 = 1 / R1;
G2 = 1 / R2;

%% Time Parameters

t_start = 0;
t_end   = 5e-3;          % 5 ms
dt      = 1e-5;          % 10 us

t = t_start:dt:t_end;

%% Preallocate Results

Vout = zeros(size(t));
Vin  = zeros(size(t));
I_Vs = zeros(size(t));

%% Initial Condition

Vout_prev = 0;

%% Transient Analysis

for k = 1:length(t)

    % Step voltage source
    if t(k) >= 0
        Vs_current = Vs;
    else
        Vs_current = 0;
    end

    % Capacitor conductance using Backward Euler
    Gc = C / dt;

    % MNA system
    %
    % Unknown vector:
    %
    % x = [Vin; Vout; I_Vs]
    %
    % The capacitor contribution includes
    % the previous time-step voltage Vout_prev.

    A = [
         G1,          -G1,             1;
        -G1,     G1 + G2 + Gc,         0;
          1,           0,              0
    ];

    b = [
        0;
        Gc * Vout_prev;
        Vs_current
    ];

    %% Solve MNA System

    x = A \ b;

    %% Store Results

    Vin(k)  = x(1);
    Vout(k) = x(2);
    I_Vs(k) = x(3);

    %% Update Previous Time-Step Voltage

    Vout_prev = Vout(k);

end

%% Display Final Result

fprintf('Transient Analysis Results\n');
fprintf('--------------------------\n');

fprintf('Input Voltage = %.2f V\n', Vs);
fprintf('R1 = %.2f Ohm\n', R1);
fprintf('R2 = %.2f Ohm\n', R2);
fprintf('C  = %.2e F\n', C);

fprintf('\nFinal Values\n');
fprintf('------------\n');

fprintf('Final Vin  = %.4f V\n', Vin(end));
fprintf('Final Vout = %.4f V\n', Vout(end));
fprintf('Final I_Vs = %.4f mA\n', I_Vs(end) * 1e3);

%% Plot Output Voltage

figure;

plot(t * 1e3, Vout, 'LineWidth', 1.5);

grid on;

xlabel('Time (ms)');
ylabel('V_{out} (V)');
title('Transient Response');

%% Plot Input and Output Voltage

figure;

plot(t * 1e3, Vin, 'LineWidth', 1.5);
hold on;
plot(t * 1e3, Vout, 'LineWidth', 1.5);

grid on;

xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Input and Output Voltage');

legend('V_{in}', 'V_{out}', 'Location', 'southeast');
