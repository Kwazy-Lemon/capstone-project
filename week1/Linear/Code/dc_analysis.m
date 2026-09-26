%% DC Analysis using Modified Nodal Analysis (MNA)
% Week 1 - Linear Circuit Prototype
%
% Circuit:
%
%        R1
% Vin ---/\/\--- Vout
%  |             |
%  |             R2
%  |             |
%  |            GND
%  |
% Voltage Source
%  |
% GND
%
% Voltage source: 10 V
% R1 = 1 kOhm
% R2 = 2 kOhm

clear;
clc;

%% Circuit Parameters

Vs = 10;           % Voltage source [V]

R1 = 1e3;          % R1 [Ohm]
R2 = 2e3;          % R2 [Ohm]

G1 = 1 / R1;
G2 = 1 / R2;

%% MNA System
%
% Unknown vector:
%
% x = [Vin; Vout; I_Vs]
%
% Vin  : voltage at input node
% Vout : voltage at output node
% I_Vs : current through the voltage source
%
% A*x = b

A = [
    G1,      -G1,          1;
    -G1,  G1 + G2,          0;
    1,       0,           0
    ];

b = [
    0;
    0;
    Vs
    ];

%% Solve MNA Equations

x = A \ b;

%% Extract Results

Vin = x(1);
Vout = x(2);
I_Vs = x(3);

%% Display Results

fprintf('DC Analysis Results\n');
fprintf('-------------------\n');

fprintf('Vin  = %.4f V\n', Vin);
fprintf('Vout = %.4f V\n', Vout);
fprintf('I_Vs = %.4f mA\n', I_Vs * 1e3);

%% Calculate Branch Currents

I_R1 = (Vin - Vout) / R1;
I_R2 = Vout / R2;

fprintf('\nCurrent Verification\n');
fprintf('--------------------\n');

fprintf('I_R1 = %.4f mA\n', I_R1 * 1e3);
fprintf('I_R2 = %.4f mA\n', I_R2 * 1e3);

fprintf('\nKCL Error = %.4e A\n', I_R1 - I_R2);
