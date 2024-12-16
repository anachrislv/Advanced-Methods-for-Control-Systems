% Parameters
N_steps = 15;
x_est_Kalman = zeros(1, N_steps); % Estimated states (EKF)
P_est = zeros(1, N_steps); % Estimated covariance
P_est(1) = 1; % Initial covariance
x_est_Kalman(1) = 0; % Initial estimate

Q = 10; % Process noise variance
R = 1;  % Measurement noise variance

% EKF Loop
for k = 1:N_steps-1
    % Prediction Step
    x_pred = 0.5 * x_est_Kalman(k) + bet * x_est_Kalman(k) / (1 + x_est_Kalman(k)^2) + 8 * cos(1.2 * k);
    F_k = 0.5 + (bet * (1 - x_est_Kalman(k)^2)) / (1 + x_est_Kalman(k)^2)^2; % Jacobian of f
    P_pred = F_k * P_est(k) * F_k' + Q;

    % Update Step
    H_k = alp + x_pred / 10; % Jacobian of h
    K_k = P_pred * H_k' / (H_k * P_pred * H_k' + R); % Kalman gain
    y_pred = alp * x_pred + x_pred^2 / 20; % Predicted measurement
    x_est_Kalman(k+1) = x_pred + K_k * (y(k) - y_pred);
    P_est(k+1) = (1 - K_k * H_k) * P_pred;
end

% Υπολογισμός MSE
MSE_Kalman = sum((x - x_est_Kalman).^2) / N_steps;

% Αποτελέσματα
fprintf('MSE for Kalman Filter: %.4f\n', MSE_Kalman);

% Plot Results
figure;
plot(1:N_steps, x, 'g', 'LineWidth', 1.5, 'DisplayName', 'True State');
hold on;
plot(1:N_steps, x_est_Kalman, 'b--', 'LineWidth', 1.5, 'DisplayName', 'EKF Estimate');
legend;
title('EKF State Estimation');
xlabel('Time Step');
ylabel('State');
grid on;
