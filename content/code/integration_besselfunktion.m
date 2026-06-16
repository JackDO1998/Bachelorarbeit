clear; close all;
%% Funktion zur berechnung der Spektralfunktion S
function S = synchrotron_S_numerisch(x)
    S = zeros(size(x));
    
    for j = 1:length(x)
        
        % Schrittweite anpassen: feiner nahe x, gröber weiter draußen
        xi = logspace(log10(x(j)), log10(100), 50000);
        % Modifizierte Besselfunktion 2. Art mit Ordnung 5/3
        K = besselk(5/3, xi);
        % Spektralfunktion des Dipols bestehend aus Vorfaktor und Integral
        % mittels Trapezregel über die modifizierte Besselfunktion K
        S(j) = (9*sqrt(3)/(8*pi)) * x(j) * trapz(xi, K);
        
    end
end

%% Eingabewerte
R = 12.2;       % Bahnradius der Elektronen im Dipolfeld 
RU = 38.2;      % Bahnradius der Elektronen im Undulator
Ib = 0.130;     % Strahlstrom
E = 1.4e9;      % Gesamtenergie Elektron


c = 299792458;  % Lichtgeschwindigkeit
E0 = 510.99895069e3;    % Ruheenergie Elektron
g = E/E0;         % Lorenzfaktor
beta = sqrt(1-(1/g^2));
e = 1.602176634e-10;    % Elementarladung
eps0 = 8.8541878188e-12;    % Elektrische Feldkonstante
lambdaU = 0.25; % Undulatorpeiode
lambda = lambdaU-lambdaU*beta*c;    % Fundamentale Undulatorlinie



x = logspace(-3, 2, 1000);
S = synchrotron_S_numerisch(x);
IntS = trapz(x,S);
fprintf('Integral von S(x) von 0,0001 bis 100 um die Normierungsbedingung Integral S(x)=1 zu prüfen %.6f\n', IntS);



% Integriere F(x) von a bis b
a = 4 * pi * RU / (3 * g^3 * (lambda-lambdaU));    % untere Grenze
b = 4 * pi * RU / (3 * g^3 * (lambda-lambdaU));     % obere Grenze

x_int = logspace(log10(a), log10(b), 5000);
S_vals = synchrotron_S_numerisch(x_int);
integral_S = trapz(x_int, S_vals);

fprintf('Integral von S(x) von %.4f bis %.4f = %.6f\n', a, b, integral_S);
Pges = Ib * e * g^4 / (3 * eps0 *R);
Pband = Pges * integral_S;


%% Plot

figure;
loglog(x, S, 'r-', 'LineWidth', 2); hold on;
xlabel('x = \omega/\omega_c');
ylabel('S(x)');
title('Spektralfunktion eines Dipols');
legend('S(x) numerisch','Location', 'northeast');
grid on;