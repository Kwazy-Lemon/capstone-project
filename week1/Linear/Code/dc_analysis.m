%% DC Analysis 
% Circuit:
%
% Vin ---R1--- Vout
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
% x = [Vin; Vout; I_Vs]
% I_Vs : current through the voltage source
% A*x = b

clear;
clc;


Vs = 10;         

R1 = 1e3;          
R2 = 2e3;          

G1 = 1 / R1;
G2 = 1 / R2;

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

x = A \ b;

Vin = x(1);
Vout = x(2);
I_Vs = x(3);


fprintf('DC Analysis Results\n');
fprintf('-------------------\n');

fprintf('Vin  = %.4f V\n', Vin);
fprintf('Vout = %.4f V\n', Vout);
fprintf('I_Vs = %.4f mA\n', I_Vs * 1e3);

I_R1 = (Vin - Vout) / R1;
I_R2 = Vout / R2;

fprintf('\nCurrent Verification\n');
fprintf('--------------------\n');

fprintf('I_R1 = %.4f mA\n', I_R1 * 1e3);
fprintf('I_R2 = %.4f mA\n', I_R2 * 1e3);

fprintf('\nKCL Error = %.4e A\n', I_R1 - I_R2);
