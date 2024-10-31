% Πίνακες συστήματος
A_ex = 1e3 * [0, -1.3333, 0; 5.7143, -1.4286, 0; 0, 0.0010, 0];
B_ex = 1e4 * [7.5; -8.0357; 0];

% Πίνακες βαρών για LQ έλεγχο
Q_ex = diag([1, 1, 1]); % Πίνακας βαρών για τις καταστάσεις
R = 1;               % Βάρος για τον έλεγχο

% Εφαρμογή της συνάρτησης icare για την επίλυση της εξίσωσης Riccati
[X, K, ~] = icare(A_ex, B_ex, Q_ex, R);

% Εμφάνιση των αποτελεσμάτων
disp('Η λύση της εξίσωσης Riccati (X):');
disp(X);
disp('Ο πίνακας ανάδρασης εισόδου (K):');
disp(K);
