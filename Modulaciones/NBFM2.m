%%%%%%%%%%%%%%%%%%%%%%% NBFM (Banda Angosta) %%%%%%%%%%%%%%%%%%%%%%

%Cargar la señal de audio
clc;clear;
rootdirectory = 'Z:\Downloads';
[ms, Fs] = audioread('AudioTarea.m4a');

m= mean(ms, 2);
m = m(:);

% Parámetros
Am = max(abs(m)); % Amplitud del mensaje
fm = 239; % Frecuencia del mensaje (Hz)
t = (0:length(m)-1) / Fs;


%Parametros de la señal portadora
fc = 1000;                  % Frecuencia de la señal portadora 1kHz
Ac = 8;                     % Amplitud de la señal portadora
c = Ac*cos(2*pi*fc*t');     % Señal portadora
c1= -(Ac*sin(2*pi*fc*t'));  % Señal portadora 90 GRADS
figure(1);
subplot(3, 2, 3);
plot(t,c);
xlabel('Tiempo (s)');
xlim([0 5/fc])
ylabel('Amplitud');
title('Señal portadora');

%Representacion de la portadora en la frecuencia
Nc = length(c);                   % Longitud de la señal c
dftc = fftshift(fft(c));          % Coloca la componente cero en el centro del espectro 
f02 = (-Nc/2:Nc/2-1)*(Fs/Nc);     % Base de frecuencias centradas en 0
DEPc = (1/(Fs*Nc))*abs(dftc).^2;  % Densidad espectral de potencia de c


subplot(3, 2, 4);
plot(f02,10*log10(DEPc))
xlim([-1.5*fc 1.5*fc]);
ylim([-80 50]);
xlabel('f [Hz]');
ylabel('|C(f)|^2 [dB]');
title('PSD de la señal portadora');
grid on

% Representación de la moduladora en el tiempo

subplot(3, 2, 1);
plot(t,m);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal moduladora');

% Representacion de la moduladora en la frecuencia
Nm = length(m);                   
dftm = fftshift(fft(m));          
f01 = (-Nm/2:Nm/2-1)*(Fs/Nm);      
DEPm = (1/(Fs*Nm))*abs(dftm).^2;  

subplot(3, 2, 2);
plot(f01,10*log10(DEPm))
xlim([-1.5*fc 1.5*fc]);
ylim([-60 20]);
xlabel('f [Hz]');
ylabel('|M(f)|^2 [dB]');
title('PSD de la señal moduladora ');
grid on

K=0.8;

%
% Señal NBFM modulada en el tiempo
nbfm_modulated = c + K .*c1.*(m - pi/2);  %
subplot(3, 2, 5);
plot(t, nbfm_modulated);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal NBFM');

Ny = length(nbfm_modulated);                   
dfty = fftshift(fft(nbfm_modulated));          
f0 = (-Ny/2:Ny/2-1)*(Fs/Ny);      
DEPy = (1/(Fs*Ny))*abs(dfty).^2;  

% PSD de la señal NBFM modulada
subplot(3, 2, 6);
plot(f0,10*log10(DEPy))
xlim([-1.5*fc 1.5*fc]);
ylim([-50 20]);
xlabel('f [Hz]');
ylabel('|S(f)|^2 [dB]');
title('PSD señal NBFM modulada')
grid on

%sound (y, Fs);
%

%{
% Señal NBFM modulada en el tiempo
nbfm_modulated = fmmod(m, fc, Fs, K);  %
subplot(3, 2, 5);
plot(t, nbfm_modulated);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal NBFM');

Ny = length(nbfm_modulated);                   
dfty = fftshift(fft(nbfm_modulated));          
f0 = (-Ny/2:Ny/2-1)*(Fs/Ny);      
DEPy = (1/(Fs*Ny))*abs(dfty).^2;  


% PSD de la señal NBFM modulada
subplot(3, 2, 6);
plot(f0,10*log10(DEPy))
xlim([-1.5*fc 1.5*fc]);
ylim([-120 0]);
xlabel('f [Hz]');
ylabel('|S(f)|^2 [dB]');
title('PSD señal NBFM')
grid on

%Demodulacion%
y = fmdemod(nbfm_modulated, fc, Fs, K);
figure (5)
subplot(3, 1, 1);
plot(t, y);
xlabel('Tiempo (s)');
ylim([-1 1]);
ylabel('Amplitud');
title('Señal Demodulada');

%}


%%%%%%%%%%%%%%%%%%%%%%%% Ruido %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Agregar ruido a las señales moduladas

snr_low = 50;       % Relación señal/ruido (bajo)
snr_medium = 15;     % Relación señal/ruido (medio)
snr_high = 2;       % Relación señal/ruido (alto)


% NBFM con ruido
nbfm_modulated_noisy_low = awgn(nbfm_modulated, snr_low);
nbfm_modulated_noisy_medium = awgn(nbfm_modulated, snr_medium);
nbfm_modulated_noisy_high = awgn(nbfm_modulated, snr_high);


figure(4);
subplot(3,1,1)
plot(t,nbfm_modulated_noisy_low);
xlim([5.56,5.565]);
ylim([-11 11]);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación NBFM con ruido bajo ')
subplot(3,1,2);
plot(t,nbfm_modulated_noisy_medium);
xlim([5.56,5.565]);
ylim([-11 11]);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación NBFM con ruido medio ')
subplot(3,1,3);
plot(t,nbfm_modulated_noisy_high);
xlim([5.56,5.565]);
ylim([-11 11]);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación NBFM con ruido alto ')
%{
%Demodulacion%
nbfm_demodulated_noisy_low1 = fmdemod(nbfm_modulated_noisy_low, fc, Fs, K);
nbfm_demodulated_noisy_medium1 = fmdemod(nbfm_modulated_noisy_medium, fc, Fs, K);
nbfm_demodulated_noisy_high1 = fmdemod(nbfm_modulated_noisy_high, fc, Fs, K);
%}


%Demodulación utilizando el discriminador de frecuencia
nbfm_demodulated_noisy_low = diff(unwrap(angle(hilbert(nbfm_modulated_noisy_low)))) * (Fs / (2 * pi * K));
nbfm_demodulated_noisy_medium = diff(unwrap(angle(hilbert(nbfm_modulated_noisy_medium)))) * (Fs / (2 * pi * K));
nbfm_demodulated_noisy_high = diff(unwrap(angle(hilbert(nbfm_modulated_noisy_high)))) * (Fs / (2 * pi * K));
nbfm_demodulated = fmdemod(nbfm_modulated,fc, Fs,K);

% Demodulación con diferentes niveles de ruido %

t_noisy_low = (0:length(nbfm_demodulated_noisy_low)-1) / Fs;
t_noisy_medium = (0:length(nbfm_demodulated_noisy_medium)-1) / Fs;
t_noisy_high = (0:length(nbfm_demodulated_noisy_high)-1) / Fs;

figure(3);
subplot(3,2,1)
plot(t_noisy_low,nbfm_demodulated_noisy_low)
title('Señal demodulada con ruido bajo')

Nb = length(nbfm_demodulated_noisy_low);                   
dftb = fftshift(fft(nbfm_demodulated_noisy_low));          
f0b = (-Nb/2:Nb/2-1)*(Fs/Nb);      
DEPb = (1/(Fs*Nb))*abs(dftb).^2;  


% PSD Señal Demodulada con ruido baj
subplot(3, 2, 2);
plot(f0b,10*log10(DEPb))
xlim([-1.5*fc 1.5*fc]);
ylim([-120 0]);
xlabel('f [Hz]');
ylabel('|S(f)|^2 [dB]');
title('PSD Señal Demodulada con ruido bajo')
grid on

subplot(3,2,3)
plot(t_noisy_medium,nbfm_demodulated_noisy_medium)
title('Señal demodulada con ruido medio')


Nm = length(nbfm_demodulated_noisy_medium);                   
dftm = fftshift(fft(nbfm_demodulated_noisy_medium));          
f0m = (-Nm/2:Nm/2-1)*(Fs/Nm);      
DEPm = (1/(Fs*Nm))*abs(dftm).^2;  


% PSD Señal Demodulada con ruido medio
subplot(3, 2, 4);
plot(f0m,10*log10(DEPb))
xlim([-1.5*fc 1.5*fc]);
ylim([-120 0]);
xlabel('f [Hz]');
ylabel('|S(f)|^2 [dB]');
title('PSD Señal Demodulada con ruido bajo')
grid on

subplot(3,2,5)
plot(t_noisy_high,nbfm_demodulated_noisy_high)
title('Señal demodulada con ruido alto')

Nh = length(nbfm_demodulated_noisy_high);                   
dfth = fftshift(fft(nbfm_demodulated_noisy_high));          
f0h = (-Nh/2:Nh/2-1)*(Fs/Nh);      
DEPh = (1/(Fs*Nh))*abs(dfth).^2;  


% PSD Señal Demodulada con ruido alto
subplot(3, 2, 6);
plot(f0h,10*log10(DEPh))
xlim([-1.5*fc 1.5*fc]);
ylim([-120 0]);
xlabel('f [Hz]');
ylabel('|S(f)|^2 [dB]');
title('PSD Señal Demodulada con ruido medio')
grid on

figure (6)
subplot(2, 1, 1);
plot(t, nbfm_demodulated );
ylim([-1000 1000 ]);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Demodulada');

Nz = length(nbfm_demodulated);
dftz = fftshift(fft(nbfm_demodulated));
f0z = (-Nz/2:Nz/2-1) * (Fs/Nz);
DEPz = (1/(Fs*Nz)) * abs(dftz).^2;

% PSD de la señal demodulada sin ruido
subplot(2,1,2)
plot(f0z,10*log10(DEPz))
xlim([-1.5*fc 1.5*fc]);
ylim([-80 20]);
xlabel('f [Hz]');ylabel('|Z(f)|^2 [dB]')
title('PSD señal demodulada sin ruido')
grid on
%}

figure (8)
subplot(3, 1, 1);
plot(t_noisy_low, nbfm_demodulated_noisy_low );
ylim([0 3000]);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Demodulada Con Ruido Bajo');

subplot(3, 1, 2);
plot(t_noisy_medium, nbfm_demodulated_noisy_medium );
ylim([0 3000]);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Demodulada Con Ruido Medio');

subplot(3, 1, 3);
plot(t_noisy_high, nbfm_demodulated_noisy_high );
ylim([0 3000]);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Demodulada Con Ruido Alto');

% Se añade este plot para comparar

% figure (9)
% subplot(2, 1, 1);
% plot(t,m);
% xlim([5.6 5.63]);
% xlabel('Tiempo (s)');
% ylabel('Amplitud');
% title('Señal moduladora');
% 
% subplot(2, 1, 2);
% plot(t, nbfm_modulated);
% xlim([5.6 5.63]);
% xlabel('Tiempo (s)');
% ylabel('Amplitud');
% title('Señal NBFM');

%soundsc(nbfm_demodulated,Fs)
%soundsc(nbfm_demodulated_noisy_low,Fs)
%soundsc(nbfm_demodulated_noisy_medium,Fs)
%soundsc(nbfm_demodulated_noisy_high,Fs)
