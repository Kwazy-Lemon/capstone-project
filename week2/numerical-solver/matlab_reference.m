function matlab_reference()
root = fileparts(mfilename('fullpath'));
out = fullfile(root, 'results_matlab');
py = fullfile(root, 'results_python');
if ~exist(out, 'dir'), mkdir(out); end

R1 = 1000; R2 = 2000; Vs = 5;
g1 = 1/R1; g2 = 1/R2;
G = sparse([g1, -g1, 1; -g1, g1+g2, 0; 1, 0, 0]);
b = [0; 0; Vs];
x = G \ b;                        % Solve G*x=b; do not use inv(G)*b.
IR1 = (x(1)-x(2))/R1;
IR2 = x(2)/R2;
valueA = [x; IR1; IR2];
exactA = [5; 10/3; -1/600; 1/600; 1/600];
resA = G*x-b;
G2 = sparse([1e-3, 0, 1; 0, 0.5e-3, -1; 1, -1, 0]);
b2 = [0; 1e-3; 2];
x2 = G2 \ b2;
valueB = [x2; x2(1)/1000; x2(2)/2000];
exactB = [4/3; -2/3; -4/3000; 4/3000; -1/3000];
resB = G2*x2-b2;

circuit = [repmat({'divider'},5,1); repmat({'floating_source'},5,1)];
quantity = {'V_in'; 'V_out'; 'I_V1_in_to_ground'; 'I_R1_in_to_out'; ...
    'I_R2_out_to_ground'; 'V_n1'; 'V_n2'; 'I_V1_n1_to_n2'; ...
    'I_R1_n1_to_ground'; 'I_R2_n2_to_ground'};
unit = {'V'; 'V'; 'A'; 'A'; 'A'; 'V'; 'V'; 'A'; 'A'; 'A'};
matlabValue = [valueA; valueB];
analyticValue = [exactA; exactB];
absError = abs(matlabValue-analyticValue);
tol = [1e-10;1e-10;1e-12;1e-12;1e-12;1e-10;1e-10;1e-12;1e-12;1e-12];
assert(all(absError <= tol), 'DC analytic check failed.');
dc = table(circuit,quantity,unit,matlabValue,analyticValue,absError);
disp(dc);
writetable(dc,fullfile(out,'dc_results_matlab.csv'));
fprintf('DC KCL residuals (A): %.3e, %.3e, %.3e, %.3e\n', ...
    resA(1),resA(2),resB(1),resB(2));
fprintf('DC source-constraint residuals (V): %.3e, %.3e\n',resA(3),resB(3));

dcPath = fullfile(py,'dc_results.csv');
if exist(dcPath,'file')
    p = readtable(dcPath);
    assert(height(p)==height(dc), 'DC row count mismatch.');
    assert(all(strcmp(p.circuit,circuit)) && all(strcmp(p.quantity,quantity)) ...
        && all(strcmp(p.unit,unit)), 'DC circuit/quantity/unit mismatch.');
    pythonValue = p.python;
    absDifference = abs(pythonValue-matlabValue);
    passed = absDifference <= tol;
    comparison = table(circuit,quantity,unit,pythonValue,matlabValue, ...
        absDifference,tol,passed);
    disp(comparison);
    writetable(comparison,fullfile(out,'dc_python_matlab_comparison.csv'));
    assert(all(passed), 'Python/MATLAB DC comparison failed.');
else
    fprintf('Python DC comparison pending: run run_benchmarks.py first.\n');
end

R = 1000; C = 1e-6; g = 1/R;
G = sparse([g,-g,1; -g,g,0; 1,0,0]);
M = sparse(3,3); M(2,2) = C;
fc = 1/(2*pi*R*C);
frequency_Hz = unique([logspace(0,5,301),fc]).';
bAC = [0;0;1];
H = complex(zeros(numel(frequency_Hz),1));
for k = 1:numel(frequency_Hz)
    omega = 2*pi*frequency_Hz(k);
    xAC = (G+1i*omega*M) \ bAC;
    H(k) = xAC(2)/xAC(1);
end
HExact = 1./(1+1i*2*pi*frequency_Hz*R*C);
assert(max(abs(H-HExact)) < 1e-12,'AC analytic check failed.');
H_real = real(H); H_imag = imag(H);
gain_dB = 20*log10(abs(H)); phase_deg = angle(H)*180/pi;
ac = table(frequency_Hz,H_real,H_imag,gain_dB,phase_deg);
writetable(ac,fullfile(out,'ac_results_matlab.csv'));

fig = figure('Name','RC AC reference');
subplot(2,1,1);
semilogx(frequency_Hz,gain_dB,'LineWidth',1.4); grid on;
ylabel('Gain (dB)'); title(sprintf('RC low-pass, fc = %.3f Hz',fc));
subplot(2,1,2);
semilogx(frequency_Hz,phase_deg,'LineWidth',1.4); grid on;
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
print(fig,fullfile(out,'ac_frequency_response_matlab.png'),'-dpng','-r160');

acPath = fullfile(py,'ac_results.csv');
if exist(acPath,'file')
    p = readtable(acPath);
    assert(height(p)==height(ac),'AC sample count mismatch.');
    assert(max(abs(p.frequency_Hz-frequency_Hz)./max(1,frequency_Hz))<1e-10, ...
        'AC frequency grids differ.');
    absComplexDifference = abs((p.H_real+1i*p.H_imag)-H);
    writetable(table(frequency_Hz,absComplexDifference), ...
        fullfile(out,'ac_python_matlab_comparison.csv'));
    assert(max(absComplexDifference)<1e-10,'Python/MATLAB AC comparison failed.');
end

Vs = 5; tau = R*C; tstop = 5*tau;
dt_over_tau = [0.2;0.1;0.05;0.025;0.0125];
dt_s = tau*dt_over_tau;
max_error_V = zeros(size(dt_s)); rmse_V = max_error_V;
steps = zeros(size(dt_s)); observed_order = nan(size(dt_s));
fig2 = figure('Name','RC transient reference');
tfine = linspace(0,tstop,1001).';
plot(tfine*1000,Vs*(1-exp(-tfine/tau)),'k','LineWidth',1.8); hold on;
labels = {'Analytical'};
for j = 1:numel(dt_s)
    dt = dt_s(j);
    steps(j) = round(tstop/dt);
    t = (0:steps(j)).'*dt;
    X = zeros(numel(t),3);
    X(1,:) = [Vs,0,-Vs/R];
    A = G+M/dt;
    for n = 2:numel(t)
        rhs = [0;0;Vs] + (M/dt)*X(n-1,:).';
        X(n,:) = (A \ rhs).';
    end
    exact = Vs*(1-exp(-t/tau));
    error = X(:,2)-exact;
    max_error_V(j) = max(abs(error));
    rmse_V(j) = sqrt(mean(error.^2));
    exactDiscrete = Vs*(1-(1+dt/tau).^(-(0:steps(j)).'));
    assert(max(abs(X(:,2)-exactDiscrete))<1e-10,'BE recurrence check failed.');
    if j > 1
        observed_order(j) = log(max_error_V(j-1)/max_error_V(j))/log(2);
    end
    if any(j==[1,3,5])
        plot(t*1000,X(:,2),'--','LineWidth',1.1);
        labels{end+1} = sprintf('BE dt = %g us',dt*1e6); %#ok<AGROW>
    end
end
grid on; xlabel('Time (ms)'); ylabel('Output voltage (V)');
title('RC step response'); legend(labels,'Location','southeast');
print(fig2,fullfile(out,'transient_comparison_matlab.png'),'-dpng','-r160');
errTable = table(dt_s,dt_over_tau,steps,max_error_V,rmse_V,observed_order);
disp(errTable);
writetable(errTable,fullfile(out,'timestep_error_matlab.csv'));
assert(all(diff(max_error_V)<0),'Smaller timesteps should reduce the RC error.');

errorPath = fullfile(py,'timestep_error.csv');
if exist(errorPath,'file')
    p = readtable(errorPath);
    assert(height(p)==height(errTable),'Timestep count mismatch.');
    assert(max(abs(p.dt_s-dt_s))<1e-14,'Timestep grids differ.');
    difference_V = abs(p.max_error_V-max_error_V);
    writetable(table(dt_s,difference_V),fullfile(out,'timestep_python_matlab_comparison.csv'));
    assert(max(difference_V)<1e-9,'Python/MATLAB timestep errors disagree.');
end
fprintf('MATLAB calculations complete. Outputs: %s\n',out);
end
