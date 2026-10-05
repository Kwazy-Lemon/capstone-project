%% RC transient analysis using MNA and Backward Euler
clear; clc; close all;

% circuit and time settings
Vs = 10; R1 = 1000; R2 = 2000; C = 100e-9;
dt = 10e-6;                       % time step: 10 us
t_end = 5e-3;                     % end time: 5 ms
t = 0:dt:t_end;
G1 = 1/R1; G2 = 1/R2; Gc = C/dt;

Vin = Vs*ones(size(t));
Vout = zeros(size(t));            % Vout(1) = 0 at t = 0
I_source = zeros(size(t));
I_source(1) = -Vs/R1;

% solve each time step, starting AFTER the initial sample
A = [G1, -G1, 1; -G1, G1+G2+Gc, 0; 1, 0, 0];
for k = 2:length(t)
    b = [0; Gc*Vout(k-1); Vs];
    x = A\b;                     % x = [Vin; Vout; source current]
    Vin(k) = x(1);
    Vout(k) = x(2);
    I_source(k) = x(3);
end

% compare with the exact RC response
Vfinal = Vs*R2/(R1+R2);
tau = (R1*R2/(R1+R2))*C;
Vexact = Vfinal*(1-exp(-t/tau));
max_error = max(abs(Vout-Vexact));

fprintf('Initial Vout = %.6f V\n',Vout(1));
fprintf('Vout at dt = %.6f V\n',Vout(2));
fprintf('Final Vout = %.6f V\n',Vout(end));
fprintf('Final source current = %.6f mA\n',I_source(end)*1e3);
fprintf('Time constant = %.6f us\n',tau*1e6);
fprintf('Maximum voltage error = %.6f V\n',max_error);

% save plots and results
folder = fullfile(fileparts(mfilename('fullpath')),'figures');
if ~exist(folder,'dir'), mkdir(folder); end

figure('Color','w');
plot(t*1e3,Vin,'b-',t*1e3,Vout,'r-','LineWidth',1.5);
grid on; xlabel('Time (ms)'); ylabel('Voltage (V)');
legend('Input','Output','Location','best');
title('Input and Output Voltage');
exportgraphics(gcf,fullfile(folder,'transient_full.png'),'Resolution',200);

figure('Color','w');
plot(t*1e6,Vout,'b-o',t*1e6,Vexact,'r--','LineWidth',1.5,'MarkerSize',3);
grid on; xlim([0 300]); xlabel('Time (us)'); ylabel('Output voltage (V)');
legend('Backward Euler','Exact solution','Location','best');
title('Transient Comparison: First 300 us');
exportgraphics(gcf,fullfile(folder,'transient_comparison.png'),'Resolution',200);

data = table(t',Vin',Vout',Vexact',I_source','VariableNames', ...
    {'Time_s','Vin_V','Vout_V','Exact_V','SourceCurrent_A'});
writetable(data,fullfile(folder,'transient_results.csv'));

