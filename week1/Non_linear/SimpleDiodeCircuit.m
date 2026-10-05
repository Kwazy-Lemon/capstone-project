
%circuit parameters
Vs = 5;         
R  = 1000;       
Is = 1e-9;     
Vt = 0.02585;   
G = 1/R;
n=1; %ideal diode

%Newton's method
maxIter = 100;
tol = 1e-9;

% unknowns: x = [v1; v2; i_vs]
% A forward-bias initial guess keeps the original Newton method but
% allows it to converge within the original 100-iteration limit.
x = [Vs; 0.4; -(Vs-0.4)/R];
error_history = zeros(maxIter,1);
for k = 1:maxIter
    v1 = x(1); v2 = x(2);

    % Diode evaluation (limit exponent to avoid overflow)
    arg = min(v2/(Vt), 40);
    Id  = Is*(exp(arg) - 1);
    gd  = Is/(n*Vt)*exp(arg);
    Ieq = Id - gd*v2;

    % Assemble MNA
    Gx = [G, -G, 1;
        -G, G+gd, 0;
        1, 0, 0];

    f = [0;
        Ieq;
        0];

    b = [0;
        0;
        Vs];

    xnew = Gx\(b-f);
    % record the same error used by the original stopping test
    error_history(k) = max(abs(xnew - x));
    
    if max(abs(xnew - x)) < tol
        x = xnew;
        fprintf('Converged in %d iterations\n', k);
        break;
    end
    x = xnew;
end

v1 = x(1); v2 = x(2); i_vs = x(3);
Id_final = Is*(exp(min(v2/(n*Vt),40)) - 1);

fprintf('v1 = %.6f V\n', v1);
fprintf('v2 = %.6f V\n', v2);
fprintf('i_vs = %.6f mA\n', i_vs*1e3);
fprintf('Diode current = %.6f mA\n', Id_final*1e3);
fprintf('Voltage across R = %.6f V\n', v1-v2);
% figure1 : working for DC operating pint.
vD_plot = linspace(0, Vs, 2000);
arg_plot = min(vD_plot/(n*Vt), 40);

Id_plot = Is*(exp(arg_plot) - 1);
IR_plot = (Vs - vD_plot)/R;

figure1 = figure('Color','w');

plot(vD_plot, Id_plot*1e3, 'r-', ...
    'LineWidth', 2.2);

hold on;

plot(vD_plot, IR_plot*1e3, 'b-', ...
    'LineWidth', 2.2);

plot(v2, Id_final*1e3, 'ko', ...
    'MarkerSize', 9, ...
    'MarkerFaceColor', 'k');

grid on;
box on;

xlabel('Diode Voltage V_D (V)');
ylabel('Current (mA)');
title('Nonlinear DC Analysis - Diode Operating Point');

legend( ...
    'Diode I-V Curve', ...
    'Resistor Load Line', ...
    sprintf('Operating Point: %.4f V, %.4f mA',v2,Id_final*1e3), ...
    'Location','best');

% Zoom into the useful forward-bias region so the intersection is clear.
xlim([0, min(Vs, max(0.8,1.5*v2))]);
ylim([0, 1.2*(Vs/R)*1e3]);

% figure2 : working for the Newton convergence 
figure2 = figure('Color','w');

semilogy(1:k, max(error_history(1:k),eps), 'mo-', ...
    'MarkerFaceColor', 'm', ...
    'MarkerEdgeColor', 'm', ...
    'LineWidth', 2.0, ...
    'MarkerSize', 6);

hold on;
yline(tol, 'g--', 'Tolerance', ...
    'LineWidth', 1.5);

grid on;
box on;

xlabel('Newton Iteration');
ylabel('Maximum Update Error');
title('Newton-Raphson Convergence');

scriptDirectory = fileparts(mfilename('fullpath'));
figureDirectory = fullfile(scriptDirectory,'diode_figures_original_style');

if ~exist(figureDirectory,'dir')
    mkdir(figureDirectory);
end

exportgraphics(figure1, ...
    fullfile(figureDirectory,'diode_dc_operating_point.png'), ...
    'Resolution',300);

exportgraphics(figure2, ...
    fullfile(figureDirectory,'newton_convergence.png'), ...
    'Resolution',300);

fprintf('\nFigures saved to:\n%s\n',figureDirectory);
