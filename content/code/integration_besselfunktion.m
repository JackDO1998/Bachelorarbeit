clear; close all; clc;

%% Physikalische Konstanten
c = 299792458;          % Lichtgeschwindigkeit [m/s]
e = 1.602176634e-19;    % Elementarladung [C]
eps0 = 8.8541878188e-12;% elektrische Feldkonstante [As/Vm]
E0 = 510.99895069e3;    % Ruheenergie eines Elektrons [eV]

%% Eingabewerte
R = 12.2;       % Dipol-Bahnradius [m]
RU = 38.2;      % Undulator-Bahnradius [m]
Ib = 0.130;     % Strahlstrom [A]
E = 1.4e9;      % Elektronenenergie [eV]
lambdaU = 0.25; % Länge der Undulatorperiode [m]
K = 3.06;       % Undulatorparameter
N_U = 5;        % Anzahl der Undulatorperioden
BPa = 90e-9;      % Untere Wellenlänge des Bandpasses      
BPb = 110e-9;      % Obere Wellenlänge des Bandpasses

%% Abgeleitete Größen
g = E / E0;     % Lorenzfaktor

% Dipol
omegac = 3 * c * g^3 / (2 * R); % kritische Kreisfrequenz des Dipols [rad/s]
lambdac = 2 * pi * c / omegac;  % kritische Wellenlänge des Dipols [m]
Pges_dipol = e * g^4 * Ib / (3 * eps0 * R); % Gesamtleistung des Dipols [W] 

% Undulator
omegac_U = 3 * c * g^3 / (2 * RU);      % kritische Kreisfrequenz des Undulators [rad/s]
lambdac_U = 2 * pi * c / omegac_U;      % kritische Wellenlänge des Undulators [m]
Pges_U = e * g^4 * Ib / (3 * eps0 * RU);% Gesamtleistung des Undulators [W]

% Resonanzwellenlängen
n_harmonics = [1, 3, 5];
lambda_n = zeros(size(n_harmonics));
for i = 1:length(n_harmonics)
    n = n_harmonics(i);
    lambda_n(i) = (lambdaU / (2 * g^2 * n)) * (1 + K^2/2);
end
x1 = lambdac_U / lambda_n(1);  % Fundamentale in x-Koordinaten

%% ============================================================
%  NORMIERTE SPEKTRALFUNKTIONEN
% ============================================================

% Dipol S(x), normiert: ∫ S dx = 1
function S = synchrotron_S_normiert(x)
    S = zeros(size(x));
    for j = 1:length(x)
        xi = logspace(log10(max(x(j),1e-6)), log10(100), 50000);
        K = besselk(5/3, xi);
        S(j) = (9*sqrt(3)/(8*pi)) * x(j) * trapz(xi, K);
        
    end
end

% Undulator S(x), normiert: ∫ S dx = 1
function S_und = undulator_S_normiert(x, x1, N_U, K, n_harmonics)
    S_und = zeros(size(x));
    
    % Relative Gewichte der Harmonischen
    P_rel = zeros(size(n_harmonics));
    for i = 1:length(n_harmonics)
        n = n_harmonics(i);
        nu1 = (n-1)/2;
        nu2 = (n+1)/2;
        z = n * K^2 / (4 + 2*K^2);
        JJ = besselj(nu1, z) - besselj(nu2, z);
        P_rel(i) = n^2 * K^2 * JJ^2;
    end
    weights = P_rel / sum(P_rel);
    
    % Linien hinzufügen
    for i = 1:length(n_harmonics)
        n = n_harmonics(i);
        x_n = n * x1;
        Delta_x = x_n / N_U;
        
        arg = (x - x_n) / Delta_x;
        sinc_vals = ones(size(arg));
        nonzero = (arg ~= 0);
        sinc_vals(nonzero) = sin(pi * arg(nonzero)) ./ (pi * arg(nonzero));
        
        S_und = S_und + weights(i) * (N_U / x_n) * sinc_vals.^2;
    end
end

%% ============================================================
%  TEST DER NORMIERUNG
% ============================================================

% Dipol
x_test = logspace(-3, 2, 10000);
S_test = synchrotron_S_normiert(x_test);
integral_dipol = trapz(x_test, S_test);
fprintf('Normierung Dipol: ∫ S(x) dx = %.6f (sollte 1 sein)\n', integral_dipol);

% Undulator
x_und_test = logspace(log10(x1*0.5), log10(5*x1*1.5), 50000);
S_und_test = undulator_S_normiert(x_und_test, x1, N_U, K, n_harmonics);
integral_und = trapz(x_und_test, S_und_test);
fprintf('Normierung Undulator: ∫ S(x) dx = %.6f (sollte 1 sein)\n', integral_und);

%% ============================================================
%  SPEKTRALE LEISTUNGSDICHTE
%  dP/dλ = P_ges * S(λ_c/λ) * (λ_c/λ²)
% ============================================================

lambda_min = min(lambdac, lambdac_U) / 100;
lambda_max = max(lambdac, lambdac_U) * 100;
lambda_vals = logspace(log10(lambda_min), log10(lambda_max), 5000);

% Dipol
x_dipol = lambdac ./ lambda_vals;
S_dipol = synchrotron_S_normiert(x_dipol);
dP_dlambda_dipol = Pges_dipol * S_dipol .* (lambdac ./ lambda_vals.^2);

% Undulator
x_und = lambdac_U ./ lambda_vals;
S_und = undulator_S_normiert(x_und, x1, N_U, K, n_harmonics);
dP_dlambda_und = Pges_U * S_und .* (lambdac_U ./ lambda_vals.^2);

% Kontrolle der Gesamtleistung
P_dipol_int = trapz(lambda_vals, dP_dlambda_dipol);
P_und_int = trapz(lambda_vals, dP_dlambda_und);
fprintf('\n--- Leistungskontrolle ---\n');
fprintf('Dipol:     ∫(dP/dλ)dλ = %.4f kW (theoretisch: %.4f kW)\n', ...
        P_dipol_int/1000, Pges_dipol/1000);
fprintf('Undulator: ∫(dP/dλ)dλ = %.4f kW (theoretisch: %.4f kW)\n', ...
        P_und_int/1000, Pges_U/1000);


%% ============================================================
%  BANDPASS-INTEGRATION (Grenzen als Wellenlängen)
%  BPa, BPb: untere und obere Wellenlänge des Bandpasses [m]
% ============================================================


% Sortieren: sicherstellen dass BPa < BPb
lambda_BP_min = min(BPa, BPb);
lambda_BP_max = max(BPa, BPb);

% Umrechnung in x = λc/λ
x_BPa_dipol = lambdac / lambda_BP_min;   % x bei kleiner Wellenlänge = großes x
x_BPb_dipol = lambdac / lambda_BP_max;   % x bei großer Wellenlänge  = kleines x

x_BPa_und = lambdac_U / lambda_BP_min;
x_BPb_und = lambdac_U / lambda_BP_max;

fprintf('\n========================================\n');
fprintf('BANDPASS-INTEGRATION\n');
fprintf('========================================\n');
fprintf('Wellenlängen-Bandpass: [%.4f nm, %.4f nm]\n', ...
        lambda_BP_min*1e9, lambda_BP_max*1e9);
fprintf('\nDipol:\n');
fprintf('  x-Bereich: [%.4f, %.4f]\n', x_BPb_dipol, x_BPa_dipol);
fprintf('Undulator:\n');
fprintf('  x-Bereich: [%.4f, %.4f]\n', x_BPb_und, x_BPa_und);

%% --- Integration über λ (physikalisch direkt) ---

% Feines λ-Gitter im Bandpass
lambda_bp = linspace(lambda_BP_min, lambda_BP_max, 10000);

% --- Dipol ---
x_dipol_bp = lambdac ./ lambda_bp;
S_dipol_bp = synchrotron_S_normiert(x_dipol_bp);
dP_dlambda_dipol_bp = Pges_dipol * S_dipol_bp .* (lambdac ./ lambda_bp.^2);

% Integral über λ → Leistung im Bandpass
P_BP_dipol = trapz(lambda_bp, dP_dlambda_dipol_bp);

% Integral über S(x) – mit Jacobi-Transformation:
% ∫ S(x) dx = ∫ S(λc/λ) * (λc/λ²) dλ
Integral_S_dipol_bp = trapz(lambda_bp, S_dipol_bp .* (lambdac ./ lambda_bp.^2));

% --- Undulator ---
x_und_bp = lambdac_U ./ lambda_bp;
S_und_bp = undulator_S_normiert(x_und_bp, x1, N_U, K, n_harmonics);
dP_dlambda_und_bp = Pges_U * S_und_bp .* (lambdac_U ./ lambda_bp.^2);

P_BP_und = trapz(lambda_bp, dP_dlambda_und_bp);
Integral_S_und_bp = trapz(lambda_bp, S_und_bp .* (lambdac_U ./ lambda_bp.^2));

%% --- Ausgabe ---

fprintf('\n--- Leistung im Bandpass ---\n');
fprintf('Dipol:     P_BP = %.4e W = %.6f kW\n', P_BP_dipol, P_BP_dipol/1000);
fprintf('           Anteil an P_ges = %.2f %%\n', 100 * P_BP_dipol / Pges_dipol);
fprintf('Undulator: P_BP = %.4e W = %.6f kW\n', P_BP_und, P_BP_und/1000);
fprintf('           Anteil an P_ges = %.2f %%\n', 100 * P_BP_und / Pges_U);

fprintf('\n--- Integral über S(x) im Bandpass ---\n');
fprintf('Dipol:     ∫ S_dipol(x) dx = %.6f\n', Integral_S_dipol_bp);
fprintf('Undulator: ∫ S_und(x)  dx = %.6f\n', Integral_S_und_bp);

% Verhältnis Undulator/Dipol im Bandpass
if P_BP_dipol > 0
    fprintf('\n--- Vergleich ---\n');
    fprintf('Verhältnis Undulator/Dipol im Bandpass: %.2f\n', P_BP_und / P_BP_dipol);
end

%% --- Plot mit Bandpass-Markierung ---

figure('Position', [100, 100, 1000, 600]);

% S(x) über λ
loglog(lambda_vals * 1e9, S_dipol, 'r-', 'LineWidth', 2); hold on;
loglog(lambda_vals * 1e9, S_und,  'b-', 'LineWidth', 1.5);

% Bandpass-Bereich schattieren
yl = ylim();
fill([lambda_BP_min*1e9, lambda_BP_max*1e9, lambda_BP_max*1e9, lambda_BP_min*1e9], ...
     [yl(1), yl(1), yl(2), yl(2)], ...
     'y', 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
     'DisplayName', sprintf('Bandpass [%.2f, %.2f] nm', lambda_BP_min*1e9, lambda_BP_max*1e9));

% Bandpass-Grenzen
xline(lambda_BP_min*1e9, 'k--', sprintf('%.2f nm', lambda_BP_min*1e9), 'LineWidth', 1.5);
xline(lambda_BP_max*1e9, 'k--', sprintf('%.2f nm', lambda_BP_max*1e9), 'LineWidth', 1.5);

% Kritische Wellenlängen
xline(lambdac * 1e9,  'r:', '\lambda_c Dipol', 'LineWidth', 1.0);
xline(lambdac_U * 1e9,'b:', '\lambda_c Undulator', 'LineWidth', 1.0);

% Harmonische
for i = 1:length(n_harmonics)
    n = n_harmonics(i);
    if n == 1
        xline(lambda_n(i)*1e9, 'g-', sprintf('\\lambda_%d', n), 'LineWidth', 1.2);
    else
        xline(lambda_n(i)*1e9, 'g:', sprintf('\\lambda_%d', n), 'LineWidth', 0.8);
    end
end

xlabel('\lambda [nm]');
ylabel('S(\lambda_c/\lambda)');
title(sprintf('Normierte Spektralfunktionen | Bandpass: %.2f – %.2f nm', ...
      lambda_BP_min*1e9, lambda_BP_max*1e9));
legend('Dipol S(x)', 'Undulator S(x)', 'Bandpass', ...
       'Location', 'southwest');
grid on;
xlim([1e-1, 1e3]);
ylim([1e-4, 1e2]);
hold off;



%% ============================================================
%  PLOT: λ gegen S(x) - Normierte Spektralfunktionen
% ============================================================

% S(x) für Dipol und Undulator (bereits aus vorherigem Abschnitt)
% x_dipol = lambdac ./ lambda_vals
% S_dipol = synchrotron_S_normiert(x_dipol)
% S_und  = undulator_S_normiert(x_und, x1, N_U, K, n_harmonics)

figure('Position', [100, 100, 1000, 600]);

loglog(lambda_vals * 1e9, S_dipol, 'r-', 'LineWidth', 2); hold on;
loglog(lambda_vals * 1e9, S_und,  'b-', 'LineWidth', 1.5);

% Kritische Wellenlängen
xline(lambdac * 1e9,  'r--', '\lambda_c Dipol', 'LineWidth', 1.2);
xline(lambdac_U * 1e9,'b--', '\lambda_c Undulator', 'LineWidth', 1.2);

% Undulator-Harmonische
for i = 1:length(n_harmonics)
    n = n_harmonics(i);
    if n == 1
        xline(lambda_n(i)*1e9, 'g-', sprintf('\\lambda_%d = %.4f nm', n, lambda_n(i)*1e9), ...
               'LineWidth', 1.5);
    else
        xline(lambda_n(i)*1e9, 'g:', sprintf('\\lambda_%d', n), ...
               'LineWidth', 1.0);
    end
end

xlabel('\lambda [nm]');
ylabel('S(\lambda_c/\lambda)');
title('Normierte Spektralfunktionen: Dipol vs. Undulator');
subtitle(sprintf('∫ S(x) dx = 1 für beide Spektren (K = %.1f, N_U = %d)', K, N_U));
legend('Dipol S(x)', 'Undulator S(x)', ...
       '\lambda_c Dipol', '\lambda_c Undulator', ...
       'Location', 'southwest');
grid on;

% y-Achse anpassen
ylim([1e-6, 1]);
hold off;

%% ============================================================
%  ZOOM-PLOT: Linear um λ₁
% ============================================================
figure('Position', [150, 150, 800, 500]);

% Bereich ±10% um λ₁
lambda_zoom = linspace(lambda_n(1)*0.9, lambda_n(1)*1.1, 2000);

% Dipol S(x) im Zoom-Bereich
x_dipol_zoom = lambdac ./ lambda_zoom;
S_dipol_zoom = synchrotron_S_normiert(x_dipol_zoom);

% Undulator S(x) im Zoom-Bereich
x_und_zoom = lambdac_U ./ lambda_zoom;
S_und_zoom = undulator_S_normiert(x_und_zoom, x1, N_U, K, n_harmonics);

plot(lambda_zoom * 1e9, S_dipol_zoom, 'r-', 'LineWidth', 2); hold on;
plot(lambda_zoom * 1e9, S_und_zoom,  'b-', 'LineWidth', 1.5);

xline(lambda_n(1) * 1e9, 'g-', sprintf('\\lambda_1 = %.4f nm', lambda_n(1)*1e9), 'LineWidth', 1.5);

xlabel('\lambda [nm]');
ylabel('S(\lambda_c/\lambda)');
title(sprintf('Zoom um \\lambda_1 (N_U = %d, K = %.1f)', N_U, K));
legend('Dipol', 'Undulator', 'Location', 'northeast');
grid on;
hold off;