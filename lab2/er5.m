% Παράμετροι συστήματος
m = 0.8; % Μάζα (kg)
g = 9.81; % Επιτάχυνση λόγω βαρύτητας (m/s^2)
Ixx = 15.67e-3; % Ροπή αδράνειας γύρω από x
Iyy = 15.67e-3; % Ροπή αδράνειας γύρω από y
Izz = 28.34e-3; % Ροπή αδράνειας γύρω από z

% Σημείο ισορροπίας
x_eq = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % Καταστάσεις ισορροπίας
u_eq = [m * g; 0; 0; 0]; % Είσοδοι ισορροπίας (U1 = m*g, άλλες ροπές = 0)

% Δυναμικό σύστημα (μη γραμμικό)
f = @(x, u) [
    x(4); % dx/dt = vx
    x(5); % dy/dt = vy
    x(6); % dz/dt = vz
    u(1)/m - g; % dvz/dt
    u(2)/Ixx; % dp/dt
    u(3)/Iyy; % dq/dt
    u(4)/Izz; % dr/dt
    zeros(5, 1) % Συμπλήρωμα για συμβατότητα
];

% Διαφορές για γραμμικοποίηση
epsilon = 1e-5; % Μικρή μεταβολή

% Υπολογισμός του A (df/dx)
n = length(x_eq); % Διαστάσεις x
A = zeros(n, n); % Μητρώο A
for i = 1:n
    x_plus = x_eq;
    x_minus = x_eq;
    x_plus(i) = x_plus(i) + epsilon;
    x_minus(i) = x_minus(i) - epsilon;
    A(:, i) = (f(x_plus, u_eq) - f(x_minus, u_eq)) / (2 * epsilon);
end

% Υπολογισμός του B (df/du)
m = length(u_eq); % Διαστάσεις u
B = zeros(n, m); % Μητρώο B
for j = 1:m
    u_plus = u_eq;
    u_minus = u_eq;
    u_plus(j) = u_plus(j) + epsilon;
    u_minus(j) = u_minus(j) - epsilon;
    B(:, j) = (f(x_eq, u_plus) - f(x_eq, u_minus)) / (2 * epsilon);
end

% Εμφάνιση αποτελεσμάτων
disp('Μητρώο A:')
disp(A)
disp('Μητρώο B:')
disp(B)
