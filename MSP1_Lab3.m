%% MSP1 Labortermin 2 - 05.11.2014

%% Aufgabe 1 - Systemidentifikation
% y(t) ist Springantowrt auf u(t) = u0 * sigma(t)
% u0 = -750 N
clc; clear all; close all;

load('bruecke_sprung.mat', '-mat');
u0 = -750;

plot(t, y);
grid on;

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

%% Aufgabe 2 - Analyse des Schwingverhaltens

% 2a) Vergleich Theorie - gemessene Werte

figure;
sys = tf(1,[mB rB kB]);
step(sys);
hold;
plot(t, y/u0, '-r');
legend('simulated', 'measured');

% 2b) Anregung durch Fußgänger
u0 = -750;
u_A = 250;
f = 1.75;
u = u0 + u_A * sin(2*pi*f*t);

hold off;
[y t] = lsim(sys, u, t);
figure;
plot(t, y);
title('step response pedestrian');
xlabel('Time (seconds)');
ylabel('Amplitude');

% Einschwingzeit:  ca. 4,5 Sekunden
% Amplitude:       0,04 m

%% 2c) Amplitudengang

[H, wout] = freqresp(sys);
H_red = squeeze(H);
Amp = abs(H_red);
plot(wout, Amp);

y_max = 1.643e-4;
y_stat = Kp;

Res_ueberh = y_max / y_stat;


%% 3 - Schwingungstilger

% 3c) Test
dT = 0.1;

mT = 25;
rT = dT; %???
%rT = (2*dT)/(omega_0);
kT = omega_0^2 * mT;

b2 = mT;
b1 = rT;
b0 = kT;
a4 = mT*mB;
a3 = (mT+mB)*rT + mT*rB;
a2 = (mT+mB)*kT + mT*kB + rT*rB;
a1 = kB*rT + kT*rB;
a0 = kT*kB;

sys_ges = tf([b2 b1 b0],[a4 a3 a2 a1 a0]);
%bode(sys_ges);

hold;

[H, wout] = freqresp(sys_ges);
H_red = squeeze(H);
Amp = abs(H_red);
plot(wout, Amp, '-r');



% 3d) eingestellt


