clc; clear; close all

K = 100;%periode echantillonnage
N = 3*K; %on limite à +/- 3K
n = (-N:N)'; 
np = (-2*N:2*N)'; %pour zero-padding

fmax = 1/K;
fc = 2*fmax; %Nyquist : fc > 2*fmax
n0 = 0.2;
dn = 0.02;
k = 500;

%% Premiere etape : Filtre d'emission h
h = sinc(2*fc*(n - 2*n)); %filtre d'emission tronque
port = cos(2*pi*n0*np);
port2 = cos(2*pi*(n0+dn)*np);
port3 = cos(2*pi*(n0-dn)*np);
delay = grpdelay(h, 1, k); %calcul du retard de groupe
delay = mean(abs(delay));

t = (0:length(h)-1)';

%% Deuxieme etape : Amplitudes a
r1 = randn(length(n), 1); %generation sequences aleatoires
r2 = randn(length(n), 1);
r3 = randn(length(n), 1);

for i = 1:K:length(n)
    a(i) = sign(r1(i));
    b(i) = sign(r2(i));
    c(i) = sign(r3(i));
end

figure(1)
stem(t, h)
grid()
title("Filtre d'émission h")
xlabel("Temps discret")
ylabel("Amplitude")

figure(2)
stem(a, 'b', 'LineWidth', 2)
hold on
stem(b, 'r')
stem(c, 'g')
grid()
title("Peigne d'impulsions a")
xlabel("Temps discret")
ylabel("Amplitude")

%% Troisieme etape : Application du filtre au signal d'amplitude
figure(3)
h = [h; zeros(length(h)-1, 1)]; %zero-padding
a = [a'; zeros(length(a)-1, 1)];
b = [b'; zeros(length(b)-1, 1)];
c = [c'; zeros(length(c)-1, 1)];

sa1 = filter(h, 1, a);

sA1 = sa1 .* port; %on multiplie a la porteuse
sA2 = sA1;
sA2(1:delay) = [];

sb1 = filter(h, 1, b);

sB1 = sb1 .* port2; %on multiplie a la porteuse
sB2 = sB1;
sB2(1:delay) = [];

sc1 = filter(h, 1, c);

sC1 = sc1 .* port3; %on multiplie a la porteuse
sC2 = sC1;
sC2(1:delay) = [];


sFinal = sA1 + sB1 + sC1;
hold on
plot(sA2, 'b') %signal module avec porteuse 
stem(a, 'b', 'LineWidth', 2)
plot(sB2, 'r')
stem(b, 'r')
plot(sC2, 'g')
stem(c, 'g')
grid()
title("Signal modulé s")
xlabel("Temps discret")
ylabel("Amplitude")
legend("Signal modulé sans retard de groupe", "peigne d'amplitude a")

figure(4)
[H, w] = freqz(h, 1, k);
tftd_h = 20*log10(abs(H));

% [S, w] = freqz(sA1, 1, k);
% tftd_s = 20*log10(abs(S));
f = w/(2*pi);

[SF, w] = freqz(sFinal, 1, k);
tftd_sf = 20*log10(abs(SF));

plot(f, tftd_h)
hold on
plot(f, tftd_sf)
grid()
title("Spectres du filtre d'émission et du signal modulé")
legend("Filtre d'émission", "Signal modulé")
xlabel("Fréquences numériques")
ylabel("Energie (dB)")

%% Quatrieme etape : Demodulation
demod = sFinal .* (port + port2 + port3);
demod2 = demod;
demod2(1:delay) = [];

[D, w] = freqz(demod, 1, k);
tftd_demod = 20*log10(abs(D));

figure(5)
plot(f, tftd_h)
hold on
plot(f, tftd_sf)
plot(f, tftd_demod)
legend("Filtre d'emission", "Signal recu", "Signal recu * porteuse")
grid()

ordre = 60;
fPB = fir1(ordre, 2*0.15, 'low', hann(ordre+1));
%ordre 27 pour dn = 0.05 et fc = 0.15
%ordre 41 pour dn = 0.04
signalFinal = filter(fPB', 1, demod);

signalFinal2 = signalFinal;
signalFinal2(1:delay) = [];

for i = 1:K:length(n)
    A(i) = sign(signalFinal2(i)); %sb car sans retard de groupe, sign pour arrondir a 1 ou -1 (ici pas besoin car parfait)
end

figure(6)
stem(A)
grid()
title("Démodulation")
xlabel("Temps discret")
ylabel("Amplitude")

%% Calcul des erreurs
seuilErreur = 10^-15;
e = abs(A-a(1:length(n))');
nbErr = 0;
for i = 1:length(e)
    if e(i) < seuilErreur
        e(i) = 0;
        nbErr = nbErr;
    else
        e(i) = e(i);
        nbErr = nbErr + 1;
    end
end

tauxErreur = nbErr/(length(e));

figure(7)
stem(e, 'x')
grid()
title("Représentation des erreurs")
xlabel("Echantillons")
ylabel("Amplitude d'erreur")