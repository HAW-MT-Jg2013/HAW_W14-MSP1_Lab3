%% MSP1 Labortermin 2 - 05.11.2014

%% Aufgabe 1 - Systemidentifikation
% y(t) ist Springantowrt auf u(t) = u0 * sigma(t)
% u0 = -750 N
clc; clear all; close all;

load('bruecke_sprung.mat', '-mat');
u0 = -750;

plot(t, y);
grid on;
title('Sprungantwort Br�cke (gemessen mit F = -750 N)');
xlabel('Zeit (Sekunden)');
ylabel('Ausschlag (m)');

t1 = 0.29;
t2 = 0.57;
t3 = 0.861;
y1 = -0.0430/u0;
y2 = -0.0116/u0;
Kp = -0.0248/u0;

Delta_1 = Kp-y1;
Delta_2 = y2-Kp;

T_e = 2*(t2-t1);
omega_e = pi/(t2-t1);

Theta = log(Delta_1/Delta_2);
d = sqrt(Theta^2 / (Theta^2+pi^2) );

omega_0 = (Theta*2) / (T_e*d);

mB = 1/(omega_0^2*Kp);
rB = 2*d/(omega_0*Kp);
kB = 1/Kp;

% Eigenfrequenz Br�cke: omega_0 = 11,28 Hz
%% Aufgabe 2 - Analyse des Schwingverhaltens
sys_B = tf(1,[mB rB kB]);

%% 2a) Vergleich Theorie - gemessene Werte
figure;
step(sys_B);
hold;
plot(t, y/u0, '-r');
legend('simulated', 'measured');

% Die Parameter des Systems passen mit ausreichender Genauigkeit.

%% 2b) Anregung durch Fu�g�nger
u0 = -750;
u_A = 250;
f = 1.75;
u = u0 + u_A * sin(2*pi*f*t);

[y_ped t] = lsim(sys_B, u, t);
figure;
plot(t, y_ped);
title('step response pedestrian');
xlabel('Time (seconds)');
ylabel('Amplitude (m)');

% Anregung des Systems mit einem Fu�g�nger:
% Einschwingzeit:  ca. 4,5 Sekunden
% Amplitude:       0,04 m

%% 2c) Amplitudengang

[H_orig, w_out] = freqresp(sys_B);
H = squeeze(H_orig);
Amp = abs(H);
plot(w_out, Amp);
title('Amplitudengang');
xlabel('Frequenz omega');
ylabel('Amplitude (m)');

y_max = 1.643e-4;
y_stat = Kp;

ResFaktor = y_max / y_stat;

% Die Resonanz�berh�hung betr�gt ca. 5.

%% 3 - Schwingungstilger

% Parameter der Br�cke siehe Aufgabe 1

%% 3c) Test des Tilgers
close all;
dT = 0.1;   % in Aufgabe gegeben
mT = 25;    % in Aufgabe gegeben

kT = mT * omega_0^2;
rT = 2*dT*kT/omega_0;

b2 = mT;
b1 = rT;
b0 = kT;
a4 = mT*mB;
a3 = (mT+mB)*rT + mT*rB;
a2 = (mT+mB)*kT + mT*kB + rT*rB;
a1 = kB*rT + kT*rB;
a0 = kT*kB;

sys_ges = tf([b2 b1 b0],[a4 a3 a2 a1 a0]);

% A3c i)
hold on;
[H_orig_ges, w_out_ges] = freqresp(sys_ges);
H_ges = squeeze(H_orig_ges);
Amp_ges = abs(H_ges);
plot(w_out_ges, Amp_ges, '-r');     % System mit Tilger
plot(w_out, Amp);                   % System ohne Tilger
legend('mit Tilger', 'ohne Tilger');
grid on;
hold off;

% A3c ii)
figure;
pzmap(sys_ges, sys_B);
legend('mit Tilger', 'ohne Tilger');
grid on;

%% A3c iii)

% Sprungantwort Br�cke mit/ ohne Tilger
figure;
[y3c_ges,t3c_ges] = step(sys_ges);  % mit Tilger
[y3c_B,t3c_B] = step(sys_B);        % ohne Tilger
plot(t3c_ges, y3c_ges*u0, t3c_B, y3c_B*u0); % skalieren und plotten
legend('mit Tilger', 'ohne Tilger');
title('Sprungantwort Br�cke (simuliert mit F = -750 N)');
xlabel('Zeit (Sekunden)');
ylabel('Amplitude (m)');
grid on;

% Bewegung Tilger relativ zur Br�cke
figure;
sys_T = tf([rT kT], [mT rT kT]);
[y3c_T,t3c_T] = lsim(sys_T, y3c_ges, t3c_ges);
y3c_Trel = y3c_ges-y3c_T;
plot(t3c_T, y3c_Trel*u0);
xlabel('Zeit (Sekunden)');
ylabel('Abstand zur Br�cke rel. zur Ruhelage (m)');
title('Sprungantwort Tilger relativ zu Br�cke (simuliert mit F = -750 N)');
grid on;

%% 3d) Tilger optimieren
% Masse mT wird beibehalten, dT und kT optimieren auf eine minimale
% Resonanz�berh�hung.
% Wir gro� ist die Eigenfrequenz des Schwingungstilgers? TODO omega_opt

mT = 25;        % beibehalten

dT_opt = 0.1;           % optimieren
kT_opt = kT; %start

omega_opt = sqrt(kT_opt / mT);
rT_opt = 2*dT_opt*kT_opt/omega_opt;

b2 = mT;
b1 = rT;
b0 = kT;
a4 = mT*mB;
a3 = (mT+mB)*rT + mT*rB;
a2 = (mT+mB)*kT + mT*kB + rT*rB;
a1 = kB*rT + kT*rB;
a0 = kT*kB;

opt_ges = tf([b2 b1 b0],[a4 a3 a2 a1 a0]);



% L�sungen:
   % 


% A3d i)
% Sprungantwort Br�cke mit/ ohne Tilger
figure;
[y3d_ges,t3d_ges] = step(opt_ges);  % mit Tilger
[y3d_B,t3d_B] = step(sys_B);        % ohne Tilger
plot(t3d_ges, y3d_ges*u0, t3d_B, y3d_B*u0); % skalieren und plotten
legend('mit Tilger', 'ohne Tilger');
xlabel('Zeit (Sekunden)');
ylabel('Amplitude (m)');
title('Sprungantwort Br�cke (simuliert mit F = -750 N)');
grid on;

% Bewegung Tilger relativ zur Br�cke
figure;
opt_T = tf([rT kT], [mT rT kT]);
[y3d_T,t3d_T] = lsim(opt_T, y3d_ges, t3d_ges);
y3d_Trel = y3d_ges-y3d_T;
plot(t3d_T, y3d_Trel*u0);
xlabel('Zeit (Sekunden)');
ylabel('Abstand zur Br�cke rel. zur Ruhelage (m)');
title('Sprungantwort Tilger relativ zu Br�cke (simuliert mit F = -750 N)');
grid on;

% A3d ii)
% TODO: Reaktion der Br�cke auf den Fu�g�nger
%       - Bewegung der Br�cke (mit/ ohne Tilger)
%       - Bewegung des Tilgers relativ zur Br�cke
%       - Wie viel Platz braucht der Tilger?



