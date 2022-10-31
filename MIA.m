clc; clear; close all

K = 100;%periode echantillonnage
N = 3*K; %on limite à +/- 3K
n = (-N:N)'; 

fmax = 1/K;
fc = 2*fmax; %Nyquist : fc > 2*fmax
k = 500;

%% Premiere etape : Filtre d'emission h
h = sinc(2*fc*(n - 2*n)); %filtre d'emission tronque
delay = grpdelay(h, 1, k); %calcul du retard de groupe
delay = mean(abs(delay));

t = (0:length(h)-1)';

%% Deuxieme etape : Amplitudes a
r = randn(length(n), 1); %generation sequences aleatoires


for i = 1:K:length(n)
    a(i) = sign(r(i));
end


figure(1)
stem(t, h)
grid()
title("Filtre d'émission h")
xlabel("Temps discret")
ylabel("Amplitude")

figure(2)
stem(a)
grid()
title("Peigne d'impulsions a")
xlabel("Temps discret")
ylabel("Amplitude")

%% Troisieme etape : Application du filtre au signal d'amplitude
figure(3)
h = [h; zeros(length(h), 1)];
a = [a'; zeros(length(a), 1)];
s = filter(h, 1, a);
sb = s;
sb(1:delay) = [];
% plot(s) %avec retard de groupe
hold on
plot(sb) %sans retard de groupe
stem(a)
grid()
title("Signal modulé s")
xlabel("Temps discret")
ylabel("Amplitude")
legend("Signal modulé sans retard de groupe", "peigne d'amplitude a")

figure(4)
[H, w] = freqz(h, 1, k);
tftd_h = 20*log10(abs(H));

[S, w] = freqz(s, 1, k);
tftd_s = 20*log10(abs(S));
f = w/(2*pi);

plot(f, tftd_h)
hold on
plot(f, tftd_s)
grid()
title("Spectres du filtre d'émission et du signal modulé")
legend("Filtre d'émission", "Signal modulé")
xlabel("Fréquences numériques")
ylabel("Energie (dB)")
%energie spectrale depend entierement filtre d'emission
%peu importe si pas gain a 0dB, car deja impossible si gain (h) et energie
%(a) constantes non = 1, ca se multiplie avec la conv
%% Quatrieme etape : Demodulation
for i = 1:K:length(n)
    A(i) = sign(sb(i)); %sb car sans retard de groupe, sign pour arrondir a 1 ou -1 (ici pas besoin car parfait)
end

figure(5)
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

figure(6)
stem(e, 'x')
grid()
title("Représentation des erreurs")
xlabel("Echantillons")
ylabel("Amplitude d'erreur")