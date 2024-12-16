% Υπολογισμός MSE
MSE_Bayesian = sum((x(1:N_steps) - x_est_Bayesian).^2) / N_steps;
MSE_PF = sum((x(1:N_steps) - x_est_PF).^2) / N_steps;
MSE_PF_noRes = sum((x(1:N_steps) - x_est_PF_noRes).^2) / N_steps;

% Εμφάνιση Αποτελεσμάτων
fprintf('MSE for Bayesian Filter: %.4f\n', MSE_Bayesian);
fprintf('MSE for Particle Filter (with resampling): %.4f\n', MSE_PF);
fprintf('MSE for Particle Filter (without resampling): %.4f\n', MSE_PF_noRes);
