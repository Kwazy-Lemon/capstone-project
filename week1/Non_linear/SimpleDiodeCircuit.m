
%circuit parameters
Vs = 5;         
R  = 1000;       
Is = 1e-9;     
Vt = 0.02585;   
G = 1/R;

%Newton's method
maxIter = 100;
tol = 1e-9;

% unknowns: x = [v1; v2; i_vs]
x = [0; 0; 0];   % initial guess

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