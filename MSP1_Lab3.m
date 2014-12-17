%% MSP1 Labortermin 3 und 4 - 14.12.2014

%% Aufgabe 1 - Systemidentifikation
% y(t) ist Sprungantwort auf u(t) = u0 * sigma(t)
% u0 = -750 N
clc; clear all; close all;

load('bruecke_sprung.mat', '-mat');
u0 = -750;

plot(t, y);
grid on;
title('Sprungantwort Bruecke (gemessen mit F = -750 N)');
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

% Eigenfrequenz der Bruecke: omega_0 = 11,28 Hz
%% Aufgabe 2 - Analyse des Schwingverhaltens
sys_B = tf(1,[mB rB kB]);

%% 2a) Vergleich Theorie - gemessene Werte
figure;
step(sys_B);
hold;
plot(t, y/u0, '-r');
legend('simuliert', 'gemessen');

% Die Parameter des Systems passen mit ausreichender Genauigkeit.

%% 2b) Anregung durch Fuﬂgaenger
u0 = -750;
u_A = 250;
f = 1.75;
u = u0 + u_A * sin(2*pi*f*t);

[y_ped t] = lsim(sys_B, u, t);
figure;
plot(t, y_ped);
title('Anregung durch Fuﬂgaenger');
xlabel('Zeit (Sekunden)');
ylabel('Amplitude (m)');

% Anregung des Systems mit einem Fuﬂgaenger:
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

% Die Resonanzueberhoehung betraegt ca. 5.

%% 3 - Schwingungstilger

% Parameter der Bruecke siehe Aufgabe 1

% 3c) Test des Tilgers
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
title('Amplitudengang');
xlabel('Frequenz omega');
ylabel('Amplitude (m)');
legend('mit Tilger', 'ohne Tilger');
grid on;
hold off;

% A3c ii)
figure;
pzmap(sys_ges, sys_B);
legend('mit Tilger', 'ohne Tilger');
grid on;

%% A3c iii)

% Sprungantwort Bruecke mit/ ohne Tilger
figure;
[y3c_ges,t3c_ges] = step(sys_ges);  % mit Tilger
[y3c_B,t3c_B] = step(sys_B);        % ohne Tilger
plot(t3c_ges, y3c_ges*u0, t3c_B, y3c_B*u0); % skalieren und plotten
legend('mit Tilger', 'ohne Tilger');
title('Sprungantwort Bruecke (simuliert mit F = -750 N)');
xlabel('Zeit (Sekunden)');
ylabel('Amplitude (m)');
grid on;

% Bewegung Tilger relativ zur Bruecke
figure;
sys_T = tf([rT kT], [mT rT kT]);
[y3c_T,t3c_T] = lsim(sys_T, y3c_ges, t3c_ges);
y3c_Trel = y3c_ges-y3c_T;
plot(t3c_T, y3c_Trel*u0);
xlabel('Zeit (Sekunden)');
ylabel('Abstand zur Bruecke rel. zur Ruhelage (m)');
title('Sprungantwort Tilger relativ zu Bruecke (simuliert mit F = -750 N)');
grid on;

%% 3d) Tilger optimieren
% Masse mT wird beibehalten, dT und kT sollen so optimiert werden, 
% dass eine minimale Resonanzueberhoehung auftritt.
% Wie gro√ü ist die Eigenfrequenz des Schwingungstilgers? 9.73 Hz

% automatisches Finden der Parameter:
% ----------------------------------

% Werte initialisieren:
Res_last = 100;
w_test = 5:0.1:25;
mT = 25;                % gegebener Wert
dT_s = 0.1;             % Startwert

% kT optimieren:
for kT_s = 2200:0.1:2500  % Bereich zum Finden der optimalen Frequenz

    omega_s = sqrt(kT_s / mT);
    rT_s = 2*dT_s*kT_s/omega_s;

    b2 = mT;
    b1 = rT_s;
    b0 = kT_s;
    a4 = mT*mB;
    a3 = (mT+mB)*rT_s + mT*rB;
    a2 = (mT+mB)*kT_s + mT*kB + rT_s*rB;
    a1 = kB*rT_s + kT_s*rB;
    a0 = kT_s*kB;

    opt_ges = tf([b2 b1 b0],[a4 a3 a2 a1 a0]);

    % Resonanzueberhoehung automatisch bestimmen:
    [H_o, w_o] = freqresp(opt_ges, w_test);
    H_o_x = squeeze(H_o);
    Amp_o = abs(H_o_x);

    y_max_o = max(Amp_o);
    y_step = step(opt_ges);
    y_stat_o = y_step(end);

    Res = y_max_o / y_stat_o;

    % Erneut durchfuehren waehrend Res sich nach wie vor vermindert
    % den letzten optimalen Wert fuer Ausgabe sichern
    if Res > Res_last
        break;
    else      
        Res_last = Res;
        kT_opt = kT_s;
        rT_opt = rT_s;
        omega_opt = omega_s;
    end
end


% Werte initialisieren:
Res_last = 100;
w_test = 5:0.1:20;

% dT optimieren
for dT_s = 0.1:0.001:1     % Bereich zum Finden der optimalen Daempfung
    
    rT_s = 2*dT_s*kT_opt/omega_opt;

    b2 = mT;
    b1 = rT_s;
    b0 = kT_opt;
    a4 = mT*mB;
    a3 = (mT+mB)*rT_s + mT*rB;
    a2 = (mT+mB)*kT_opt + mT*kB + rT_s*rB;
    a1 = kB*rT_s + kT_opt*rB;
    a0 = kT_opt*kB;

    opt_ges = tf([b2 b1 b0],[a4 a3 a2 a1 a0]);

    % Resonanzueberhoehung automatisch bestimmen:
    [H_o, w_o] = freqresp(opt_ges, w_test);
    H_o_x = squeeze(H_o);
    Amp_o = abs(H_o_x);

    y_max_o = max(Amp_o);
    y_step = step(opt_ges);
    y_stat_o = y_step(end);

    Res = y_max_o / y_stat_o;

    % Erneut durchfuehren waehrend Res sich nach wie vor vermindert
    % den letzten optimalen Wert fuer Ausgabe sichern
    if Res > (Res_last+0.01)
        break;
    else      
        Res_last = Res;
        rT_opt = rT_s;
        dT_opt = omega_opt*rT_opt/(kT_opt*2);
    end
     
end

% System neuberechnen mit optimierten Werten:
% ---------------------------------------

% fuer manuelle Selektion:
% dT_s = 0.188;
% rT_opt = 2*dT_s*kT_opt/omega_opt;

b2 = mT;
b1 = rT_opt;
b0 = kT_opt;
a4 = mT*mB;
a3 = (mT+mB)*rT_opt + mT*rB;
a2 = (mT+mB)*kT_opt + mT*kB + rT_opt*rB;
a1 = kB*rT_opt + kT_opt*rB;
a0 = kT_opt*kB;

opt_ges = tf([b2 b1 b0],[a4 a3 a2 a1 a0]);

% Resonanzueberhoehung automatisch bestimmen:
[H_o, w_o] = freqresp(opt_ges);
H_o_x = squeeze(H_o);
Amp_o = abs(H_o_x);

y_max_o = max(Amp_o);
y_step = step(opt_ges);
y_stat_o = y_step(end);

Res = y_max_o / y_stat_o;

plot(w_o, Amp_o);
title('Amplitudengang');
xlabel('Frequenz omega');
ylabel('Amplitude (m)');
ylim([0 1.4e-4]);
grid on;


% Loesungen:
% omega_opt = 9.7334
% kT_opt = 2.3685e+03
% rT_opt = 87.6010
% dT_opt = 0.1800

%% A3d i)
% Sprungantwort Bruecke mit/ ohne Tilger
figure;
[y3d_ges,t3d_ges] = step(opt_ges);  % mit Tilger
[y3d_B,t3d_B] = step(sys_B);        % ohne Tilger
plot(t3d_ges, y3d_ges*u0, t3d_B, y3d_B*u0); % skalieren und plotten
legend('mit Tilger', 'ohne Tilger');
xlabel('Zeit (Sekunden)');
ylabel('Amplitude (m)');
title('Sprungantwort Bruecke (simuliert mit F = -750 N)');
grid on;

% Bewegung Tilger relativ zur Bruecke
figure;
opt_T = tf([rT_opt kT_opt], [mT rT_opt kT_opt]);
[y3d_T,t3d_T] = lsim(opt_T, y3d_ges, t3d_ges);
y3d_Trel = y3d_ges-y3d_T;
plot(t3d_T, y3d_Trel*u0);
xlabel('Zeit (Sekunden)');
ylabel('Abstand zur Bruecke rel. zur Ruhelage (m)');
title('Sprungantwort Tilger relativ zu Bruecke (simuliert mit F = -750 N)');
grid on;

%% A3d ii)
% Reaktion der Bruecke auf den Fuﬂgaenger
% - Bewegung der Bruecke (mit/ ohne Tilger)
% - Bewegung des Tilgers relativ zur Bruecke
% - Wie viel Platz braucht der Tilger?

% Anregung durch Fuﬂgaenger:
u0 = -750;
u_A = 250;
f = 1.75;
u = u0 + u_A * sin(2*pi*f*t);

% Amplitude
figure;
[y_ped t] = lsim(sys_B, u, t);
[y_ped_opt t] = lsim(opt_ges, u, t);
plot(t, y_ped, t, y_ped_opt);
title('Anregung durch Fuﬂgaenger');
xlabel('Time (seconds)');
ylabel('Amplitude (m)');
legend('ohne Tilger', 'mit optimiertem Tilger');

% Tilger relativ zu Bruecke
figure;
opt_T = tf([rT_opt kT_opt], [mT rT_opt kT_opt]);  % optimierter Tilger
[y_ped_opt t] = lsim(opt_ges, u, t);              % Bruecke+Tilger
[y3d_T t3d_T] = lsim(opt_T, y_ped_opt, t);  % Tilger: u := Bruecke

y3d_Trel = y_ped_opt-y3d_T;
plot(t3d_T, y3d_Trel);
xlabel('Zeit (Sekunden)');
ylabel('Abstand zur Bruecke rel. zur Ruhelage (m)');
title('Bewegung Tilger relativ zu Bruecke (simuliert mit Fuﬂgaenger)');
grid on;


% Platz fuer Tilger
t_min = min(y3d_Trel);
t_max = max(y3d_Trel);

% Der Tilger benoetigt einen Mindestabstand von ca. 5,5 cm zur Bruecke.

