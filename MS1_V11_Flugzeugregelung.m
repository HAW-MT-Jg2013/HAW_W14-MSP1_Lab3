%% Sehr einfache Regelung der Rollbewegung eines Flugzeugs

% Modell des Flugzeugs: alpha=2, gamma = 0.1
alpha = 1/2; 
gamma = 1;
fl = tf(alpha, [1 gamma 0]);
subplot(2,1,1);
step(fl);
subplot(2,1,2); 
impulse(fl);



%% Regelkreis mit einfachem Proportionalregler 
% Die Strecke hat I-Verhalten, daher sollte ein P-Regler gut funktionieren.
% Kp sollte man variieren, um zu sehen, wie sich das 
% Verhalten des Regelkreises dem Regelparameter aendert

% close all;

Kp = 1;     
P_Regler = tf(Kp);    % Der P-Regler
Regelkreis = feedback(P_Regler * fl, 1);
subplot(2,1,1);
step(Regelkreis);
hold on; 
impulse(Regelkreis); 
hold off;
legend('Sprungantwort', 'Impulsantwort'); 
title('Sprung- und Impulsantwort des Regelkreises'); 
subplot(2,1,2);
pzmap(Regelkreis);
ylim([-1.5, 1.5]);
xlim([-1, 0]);
title('Polstellen des Regelkreises'); 


%% Beipiel fuer Reglerauslegung 
% Finde Kp so, dass die Anstiegszeit hoechstens 1.5s ist. 
for Kp = 0.1:0.1:10 
   P_Regler = tf(Kp);    
   Regelkreis = feedback(P_Regler * fl, 1);
   [y,t] = step(Regelkreis); 
   s = stepinfo(y,t);
   if s.RiseTime < 1.5
      break;
   end
end
ltiview(Regelkreis);


%% Wurzelortskurve
% Bei Mausklick auf die Kurve werden weitere Punkte angezeigt

rlocus(fl); 
ylim([-1, 1]);



%% Simulation mit Rechtecksignal am Eingang
Kp = 6.24;     
P_Regler = tf(Kp);    
Regelkreis = feedback(P_Regler * fl, 1);

T = 40; 
t = linspace(0,T);
u = (t>T/4) & (t<3*T/4); 
yr = lsim(Regelkreis, u, t);
subplot(2,1,1); 
plot(t, yr, t, u);
legend('Ist', 'Soll'); 
xlabel('Zeit t'), ylabel('Rollwinkel'); 
title('Soll- und Istwerte'); 
% Berechnung des Querruder-Ausschlages
subplot(2,1,2); 
yq = lsim(feedback(P_Regler, fl), u, t); 
plot(t,yq); 
xlabel('Zeit t'), ylabel('Querruder'); 
title('Stellgroesse'); 



























%% Plot fuers Skript
subplot(1,2,1);
h = rlocusplot(fl);
p = getoptions(h); % get options for plot
p.Title.String = 'Wurzelortskurve'; % change title in options
p.Title.FontSize = 12;
p.XLabel.String = 'Re(s)';
p.XLabel.FontSize = 10;
p.YLabel.String = 'Im(s)';
p.YLabel.FontSize = 10;
setoptions(h,p); % apply options to plot  

c = get(gca, 'Children')
for d=c 
   cc = get(d, 'Children')
   for dd = cc{1}
      if length(dd) > 0
         set(dd,  'LineWidth', 1.5);
         
      end
   end
end


subplot(1,2,2); 
Kp = 5;     
P_Regler = tf(Kp);    % Der R-Regler
Regelkreis = feedback(P_Regler * fl, 1);
t = linspace(0,50,1000);
[y1,t] = step(Regelkreis,t);
Kp = 1;     
P_Regler = tf(Kp);    % Der R-Regler
Regelkreis = feedback(P_Regler * fl, 1);
t = linspace(0,50,1000);
[y2,t] = step(Regelkreis,t);
Kp = 0.4;     
P_Regler = tf(Kp);    % Der R-Regler
Regelkreis = feedback(P_Regler * fl, 1);
t = linspace(0,50,1000);
[y3,t] = step(Regelkreis,t);

plot(t,y1, t,y2, t,y3,  'LineWidth', 1.5);
xlim([0,50]);
xlabel('Zeit t (Sek.)'); 
ylabel('Rollwinkel');
title('Sprungantwort des Regelkreises', 'FontSize', 12);
legend('K_p = 5', 'K_p=1', 'K_p=0.4');
