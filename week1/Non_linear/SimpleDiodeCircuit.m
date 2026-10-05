%% resistor-diode DC analysis using MNA and Newton's method
clear; clc; close all;

% circuit and solver settings
Vs = 5; R = 1000; Is = 1e-9; Vt = 0.02585; n = 1;
G = 1/R;
maxIter = 100;
voltage_tol = 1e-9;
current_tol = 1e-12;
vD_initial = 0.4;

% x = [source voltage; diode voltage; source current]
x = [Vs; vD_initial; -(Vs-vD_initial)/R];
update_tol = [voltage_tol; voltage_tol; current_tol];
residual_tol = [current_tol; current_tol; voltage_tol];

% history columns: scaled update, scaled residual, KCL, voltage residual
history = zeros(maxIter,4);
iterations = 0;
converged = false;
status = 'Maximum iterations reached';

% Newton iterations
for k = 1:maxIter
    if any(~isfinite(x)) || x(2)/(n*Vt) > 700
        status = 'Unsafe state';
        break;
    end

    % replace the diode with its linear model at the current voltage
    Id = Is*expm1(x(2)/(n*Vt));
    gd = Is/(n*Vt)*exp(x(2)/(n*Vt));
    Ieq = Id - gd*x(2);
    A = [G, -G, 1; -G, G+gd, 0; 1, 0, 0];
    b = [0; -Ieq; Vs];
    step = (A\b) - x;
    if any(~isfinite(step))
        status = 'Invalid Newton step';
        break;
    end

    % limit the diode-voltage change to 0.1 V per iteration
    scale = min(1, 0.1/max(abs(step(2)),eps));
    xnew = x + scale*step;
    if any(~isfinite(xnew)) || xnew(2)/(n*Vt) > 700
        status = 'Unsafe trial state';
        break;
    end

    % check the original circuit equations at the new solution
    Id = Is*expm1(xnew(2)/(n*Vt));
    residual = [G*(xnew(1)-xnew(2)) + xnew(3); ...
                G*(xnew(2)-xnew(1)) + Id; ...
                xnew(1)-Vs];
    if any(~isfinite(residual))
        status = 'Invalid residual';
        break;
    end
    update_error = max(abs(xnew-x)./update_tol);
    residual_error = max(abs(residual)./residual_tol);
    history(k,:) = [update_error, residual_error, ...
                   max(abs(residual(1:2))), abs(residual(3))];
    x = xnew;
    iterations = k;
    if update_error <= 1 && residual_error <= 1
        converged = true;
        status = 'Converged';
        break;
    end
end

% print results (a failed run only reports its last iterate)
fprintf('Status: %s | Iterations: %d\n',status,iterations);
fprintf('Last iterate: v1 = %.9f V, vD = %.9f V, I_source = %.9f mA\n',x(1),x(2),x(3)*1e3);
if ~converged
    warning('Diode:NotConverged','No validated operating point: %s',status);
else
    fprintf('Diode current = %.9f mA\n',Id*1e3);
end

% save diagnostics and plots
folder = fullfile(fileparts(mfilename('fullpath')),'figures');
if ~exist(folder,'dir'), mkdir(folder); end
if iterations > 0
    history = history(1:iterations,:);
    fprintf('Scaled update = %.3e | Scaled residual = %.3e (both must be <= 1)\n',history(end,1:2));
    fprintf('KCL residual = %.3e A | Voltage residual = %.3e V\n',history(end,3:4));
    figure('Color','w');
    semilogy(1:iterations,max(history(:,1:2),eps),'-o','LineWidth',1.5);
    hold on; yline(1,'k--'); grid on;
    xlabel('Iteration'); ylabel('Error divided by tolerance');
    legend('Scaled update','Scaled nonlinear residual','Threshold','Location','best');
    title(['Newton Diagnostics: ',status]);
    exportgraphics(gcf,fullfile(folder,'diode_convergence.png'),'Resolution',200);
    data = array2table([(1:iterations)',history],'VariableNames', ...
        {'Iteration','ScaledUpdate','ScaledResidual','KCL_A','VoltageResidual_V'});
    writetable(data,fullfile(folder,'diode_diagnostics.csv'));
end

if converged
    vd = linspace(min(0,x(2)-0.2),max(0.8,x(2)+0.1),1000);
    figure('Color','w');
    plot(vd,Is*expm1(vd/(n*Vt))*1e3,'r-',vd,(Vs-vd)/R*1e3,'b-','LineWidth',1.5);
    hold on; plot(x(2),Id*1e3,'ko','MarkerFaceColor','k'); grid on;
    ylim([0,1.2*max([abs(Vs/R),abs(Id),eps])*1e3]);
    xlabel('Diode voltage (V)'); ylabel('Current (mA)');
    legend('Diode I-V','Resistor load line','Operating point','Location','best');
    title('Validated Diode DC Operating Point');
    exportgraphics(gcf,fullfile(folder,'diode_operating_point.png'),'Resolution',200);
end
