%% AC Analysis
% Circuit:
%
%                 R1 = 1 kOhm
% Vin  ---------R1--------o Vout
%                         |
%                         R2---- R2 = 2 kOhm ---- GND
%                         |
%                         C---- C = 100 nF ----- GND
%
% AC voltage source: 10 V
% The circuit is analyzed over a range of frequencies.
% x = [Vin; Vout; I_Vs]
%
% A*x = b

clear;
clc;
close all;

Vs = 10;               

R1 = 1e3;                
R2 = 2e3;               
C  = 100e-9;             

G1 = 1 / R1;
G2 = 1 / R2;

f = logspace(0, 6, 1000);    % 1 Hz to 1 MHz
w = 2 * pi * f;

Vout = zeros(size(f));
I_Vs = zeros(size(f));


for k = 1:length(f)

    Yc = 1j * w(k) * C;
    A = [
         G1,       -G1,             1;
        -G1,   G1 + G2 + Yc,        0;
          1,        0,              0
    ];

    b = [
        0;
        0;
        Vs
    ];

    x = A \ b;

    Vout(k) = x(2);
    I_Vs(k) = x(3);

end

Vout_mag = abs(Vout);
Vout_phase = angle(Vout) * 180 / pi;

figure;

semilogx(f, Vout_mag, 'LineWidth', 1.5);

grid on;

xlabel('Frequency (Hz)');
ylabel('|V_{out}| (V)');
title('AC Frequency Response - Magnitude');


figure;

semilogx(f, Vout_phase, 'LineWidth', 1.5);

grid on;

xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('AC Frequency Response - Phase');

fprintf('AC Analysis Results\n');
fprintf('-------------------\n');

fprintf('Input Voltage = %.2f V\n', Vs);
fprintf('R1 = %.2f Ohm\n', R1);
fprintf('R2 = %.2f Ohm\n', R2);
fprintf('C  = %.2e F\n', C);

fprintf('\nSelected Frequency Points\n');
fprintf('-------------------------\n');

selected_frequencies = [1, 10, 100, 1e3, 1e4, 1e5, 1e6];

for k = 1:length(selected_frequencies)

    [~, idx] = min(abs(f - selected_frequencies(k)));

    fprintf( ...
        'f = %8.1f Hz : |Vout| = %8.4f V, Phase = %8.2f deg\n', ...
        f(idx), ...
        Vout_mag(idx), ...
        Vout_phase(idx) ...
    );

end
