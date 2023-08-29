syms Kc % variabile simbolica da trovare nell'equazione con vpasolve
h=2 % numero di integratori dato l'ingresso u(t)
Kd=3 % ce lo da il prof
emax=0.003334 % e<=0.003334
Kp=18 % cooefficiente di P(s)
eq1= Kd^2/(Kc*Kp)*2==emax % sempre la stessa tranne il due che è il coefficiente di u(t)
Kc=double(vpasolve(eq1,Kc)) % valore associato alla variabile simbolica trovato nell'equazione 
Kc=300 % approssimazione

s=tf([1 0],1) % s
P=Kp*(s^2/40^2+1.8/40*s+1)*(s/400+1)/((s/20+1)*(s/200+1)^2) % P(s), Kp lo leggo dal testo di P(s)
C=1/s^2 % s^-h
rlocus(C*P/Kd) % luogo pos
rlocus(-C*P/Kd) % luogo neg

finestra=20*log10(Kd/0.6) % 20 * log10(x) = x|dB ; 0.6 (emax) e kd sono dati dal testo
wmax=20 % dato dal testo

%figure
wmin=1; % lo dobbiamo mettere noi, di solito sempre 1
margin(P*C*Kc/Kd,{wmin,10000}),grid on  % lo dobbiamo mettere noi
hold on  % lo dobbiamo mettere noi
a=gcf;
a.CurrentAxes=a.Children(3); % per disegnare nel subplot giusto
% parametri: x iniziale, y iniziale, ampiezza, altezza
r=rectangle('Position',[wmin,0,wmax-wmin,finestra]);  % i parametri li dobbiamo mettere noi
r.FaceColor='red';
hold off  % lo dobbiamo mettere noi

m=8; % m > 1
tau=6/100 % omega della carta (quella relativa alla rete corretrice) / omega di bode (quella che fornisce il prof)
R=(s*tau+1)/(s*tau/m+1) % m a denominatore di R(s) ANTICIPATRICE, altrimenti attenuatrice

% copia e incolla di quello di prima, aggiungo Bode a ciclo aperto all'inizio e Nichols a
% alla fine

wmin=1;
bode(P*C*Kc/Kd,{wmin,10000}),grid on
hold on
margin(R*P*C*Kc/Kd,{wmin,10000}),grid on % qui c'è R(s)


a=gcf;
a.CurrentAxes=a.Children(3); % per disegnare nel subplot giusto
% parametri: x iniziale, y iniziale, ampiezza, altezza
r=rectangle('Position',[wmin,0,wmax-wmin,finestra]);
r.FaceColor='red';
hold off

nichols(R*P*C*Kc/Kd),ngrid % ngrid -> nichols ; qui c'è R(s)

W=feedback(R*P*C*Kc,1/Kd) % catena diretta e di controreazione
step(W) % W(s) con gradino
stepinfo(W) % info W(s) con gradino
bode(W) % Bode di W(s)

% margine di fase e omega di taglio: lo leggiamo dal grafico superiore dato da margin(R*P*C*Kc/Kd,{wmin,10000}),grid on.  
% modulo alla risonanza: differenza valore di picco e valore iniziale.
% overshoot: su stepinfo.
% risetime: su stepinfo.
% banda passante: ascissa corrispondente al valore di regime in dB -3 dB (Bode).

syms d % var simbolica come prima...
wt=99.9 % 100
vpasolve(wt*d==(64-10)*pi/180) % 64 gradi è il margine di fase, 10 lo da il testo e w di taglio con il margine di fase lo leggiamo immediatamente sopra; * pi / 180 per convertire in radianti

% ZOH con Shannon
% N.B. = x minimo 10 solitamente e se è un x buono ricostruiremo bene il
% segnale con ingresso a gradino (grafico step response).
w_c1=2*178 % doppio della banda passante
Tc1=2*pi/w_c1 % Tc = 2pi/Wc
zoh1=(1-exp(-s*Tc1))/s; % ZOH(S) = (1 - e^(-s * Tc1)) / s

w_c2=16*178 % lo scriviamo dopo, 16 è il valore giusto per regolare ; x = 16
Tc2=2*pi/w_c2 % ricalcolo Tc
zoh2=(1-exp(-s*Tc2))/s; % ricalcolo ZOH(s)

% confronto tramite Bode i 2 risultati :

bode(zoh1)
hold on
%---------
bode(zoh2)
hold off

% copia e incolla di quello sopra, ma aggiungo * zoh1/Tc1 con margin

% NON necessario
bode(R*P*C*Kc/Kd,{wmin,10000})
hold on
margin(R*P*C*Kc/Kd*zoh1/Tc1,{wmin,10000}),grid on

a=gcf;
a.CurrentAxes=a.Children(3); % per disegnare nel subplot giusto
% parametri: x iniziale, y iniziale, ampiezza, altezza
r=rectangle('Position',[wmin,0,wmax-wmin,finestra]);
r.FaceColor='red';
hold off

 
% copia e incolla di quello sopra, ma aggiungo * zoh2/Tc2 con margin

bode(R*P*C*Kc/Kd,{wmin,10000})
hold on
margin(R*P*C*Kc/Kd*zoh2/Tc2,{wmin,10000}),grid on

a=gcf;
a.CurrentAxes=a.Children(3); % per disegnare nel subplot giusto
% parametri: x iniziale, y iniziale, ampiezza, altezza
r=rectangle('Position',[wmin,0,wmax-wmin,finestra]);
r.FaceColor='red';
hold off

% W(s) aggiungendo *zoh1/Tc1 e *zoh2/Tc2 e poi li confronto con ingresso a
% gradino

W1=feedback(R*P*C*Kc*zoh1/Tc1,1/Kd)
W2=feedback(R*P*C*Kc*zoh2/Tc2,1/Kd)

step(W)
hold on
step(W1)
step(W2)
hold off
legend
% N.B = guardare la richiesta e il pshase margin di ZOH1(S) e ZOH2(S)
% N.B. = step nel tempo

% Discretizzazione controllore al variare di Tc

Cs = R * Kc * C % controllore, guarda il sistema ^
Cz1 = c2d(Cs, Tc1, 'tustin')
bode(Cz1, 'g')
hold on
bode(Cs, 'b')
hold off
