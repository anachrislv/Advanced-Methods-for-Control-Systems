% Boost Converter Equilibrium Point Calculation

% Parameters
Vin_val = 9;          % Input Voltage (V)
L_val = 0.3e-3;       % Inductance (H)
C_val = 0.07e-3;      % Capacitance (F)
R_val = 10;           % Load Resistance (Ohm)
D_range = 0.1:0.1:0.9; % Range of duty cycle D

% Preallocate arrays for equilibrium points
iL_eq_vals = zeros(size(D_range));
Vc_eq_vals = zeros(size(D_range));

% Loop over each duty cycle value
for idx = 1:length(D_range)
    D_val = D_range(idx);
    
    % Solve for equilibrium conditions
    % Using the equations:
    % iL_eq = Vo / (R * (1 - D))
    % Vo_eq = Vin / (1 - D)
    
    Vc_eq_vals(idx) = Vin_val / (1 - D_val);         % Equilibrium capacitor voltage
    iL_eq_vals(idx) = Vc_eq_vals(idx) / (R_val * (1 - D_val)); % Equilibrium inductor current
end

% Display results
disp('Duty Cycle (D)   iL_eq_vals (A)    Vc_eq_vals (V)')
for idx = 1:length(D_range)
    fprintf('     %.2f        %.4f       %.4f\n', D_range(idx), iL_eq_vals(idx), Vc_eq_vals(idx));
end

% Plot the results
figure;
subplot(2,1,1);
plot(D_range, iL_eq_vals, '-o');
xlabel('Duty Cycle (D)');
ylabel('Equilibrium Current iL_{eq} (A)');
title('Equilibrium Inductor Current vs Duty Cycle');

subplot(2,1,2);
plot(D_range, Vc_eq_vals, '-o');
xlabel('Duty Cycle (D)');
ylabel('Equilibrium Voltage Vc_{eq} (V)');
title('Equilibrium Capacitor Voltage vs Duty Cycle');
