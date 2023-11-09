%%%%%%%%%%%%%%%%%%%%%%% WBFM %%%%%%%%%%%%%%%%%%%%%%

%Cargar la señal de audio
clc; clear; close;
[ms, Fs] = audioread('AudioTarea.m4a');


m = mean(ms, 2); % Señal monoaural
m = m(:);

% Parámetros

t = (0:length(m)-1) / Fs;
Am = max(abs(m)); % Amplitud máxima de la señal de audio.
fm = 239; %PARA ESTE AUDIO, frecuencia de la moduladora. 

%Parametros de la señal portadora
fc = 1000;                  % Frecuencia de la señal portadora 1kHz
Ac = 10;                     % Amplitud de la señal portadora
c = Ac*cos(2*pi*fc*t');     % Señal portadora
c1= -(Ac*sin(2*pi*fc*t'));  % Señal portadora 90 GRADS
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
xlim([-2*fc 2*fc]);
ylim([-80 50]);
xlabel('f [Hz]');
ylabel('|C(f)|^2 [dB]');
title('Señal portadora en la frecuencia ');
grid on

% Representación de la moduladora en el tiempo
figure(1);
subplot(3, 2, 1);
plot(t,m);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal moduladora');

%Representacion de la moduladora en la frecuencia
Nm = length(m);                   
dftm = fftshift(fft(m));          
f01 = (-Nm/2:Nm/2-1)*(Fs/Nm);      
DEPm = (1/(Fs*Nm))*abs(dftm).^2;  

subplot(3, 2, 2);
plot(f01,10*log10(DEPm))
xlim([-2*fc 2*fc]);
ylim([-130 -40]);
xlabel('f [Hz]');
ylabel('|M(f)|^2 [dB]');
title('Señal moduladora en la frecuencia ');
grid on

%Modulación
beta = 10; %Índice de modulación

% Señal NBFM modulada en el tiempo
s = Ac + beta .*c1.*(m - pi/2);
subplot(3, 2, 5);
plot(t, s);
ylim([-100 100]);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal WBFM');

Ny = length(s);                   
dfty = fftshift(fft(s));          
f0 = (-Ny/2:Ny/2-1)*(Fs/Ny);      
DEPy = (1/(Fs*Ny))*abs(dfty).^2;  

% PSD de la señal NBFM modulada
subplot(3, 2, 6);
plot(f0,10*log10(DEPy))
xlim([-2*fc 2*fc]);
ylim([-100 50]);
xlabel('f [Hz]');
ylabel('|S(f)|^2 [dB]');
title('PSD señal WBFM')
grid on

%Demodulación
% Señal demodulada
envelop = abs(hilbert(s)); % Envolvente de la señal modulada
lowPassFilter  = designfilt('lowpassfir', 'FilterOrder', 100, 'CutoffFrequency', 2*fm/Fs, 'SampleRate', Fs); % Filtro pasa bajos
s_demodulada = filter(lowPassFilter, envelop); % Señal demodulada
s_normalizado = s_demodulada / max(abs(s_demodulada)); % Señal normalizada

figure (6)
subplot(2, 1, 1);
plot(t, s_normalizado);
ylim([0.7 1]);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Demodulada');

Nz = length(s_normalizado);
dftz = fftshift(fft(s_normalizado));
f0z = (-Nz/2:Nz/2-1) * (Fs/Nz);
DEPz = (1/(Fs*Nz)) * abs(dftz).^2;


% PSD de la señal demodulada sin ruido
subplot(2,1,2)
plot(f0z,10*log10(DEPz))
xlim([-1.5*fc 1.5*fc]);
ylim([-120 -50]);
xlabel('f [Hz]');ylabel('|Z(f)|^2 [dB]')
title('PSD señal demodulada sin ruido')
grid on


%%%%%%%%%%%%%%%%%%%%%%%% Ruido %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nivel_ruido_bajo = 0.01;
nivel_ruido_medio = 0.1;
nivel_ruido_alto = 0.5;

ruido_bajo = nivel_ruido_bajo*randn(size(m));
ruido_medio = nivel_ruido_medio*randn(size(m));
ruido_alto = nivel_ruido_alto*randn(size(m));

audio_noisy_low = m + ruido_bajo;
audio_noisy_medium = m + ruido_medio;
audio_noisy_high = m + ruido_alto;


figure(2);
subplot(3,1,1)
plot(t,audio_noisy_low);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Mensaje con ruido bajo ')
subplot(3,1,2);
plot(t,audio_noisy_medium);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Mensaje con ruido medio ')
subplot(3,1,3);
plot(t,audio_noisy_high);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Mensaje con ruido alto ')

% audio con ruido modulado
audio_noisy_low_modulated = Ac + beta .*c1.*(audio_noisy_low - pi/2);
audio_noisy_medium_modulated = Ac + beta .*c1.*(audio_noisy_medium - pi/2);
audio_noisy_high_modulated = Ac + beta .*c1.*(audio_noisy_high - pi/2);

figure(4);
subplot(3,1,1)
plot(t,audio_noisy_low_modulated);
xlim([5 5.005]);
ylim([-200 200]);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación WBFM con ruido bajo ')
subplot(3,1,2);
plot(t,audio_noisy_medium_modulated);
xlim([5 5.005]);
ylim([-200 200]);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación WBFM con ruido medio ')
subplot(3,1,3);
plot(t,audio_noisy_high_modulated);
xlim([5 5.005]);
ylim([-200 200]);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación WBFM con ruido alto ')


% Señal demodulada con ruido
envelop_low = abs(hilbert(audio_noisy_low_modulated)); % Envolvente de la señal modulada
envelop_medium = abs(hilbert(audio_noisy_medium_modulated)); % Envolvente de la señal modulada
envelop_high = abs(hilbert(audio_noisy_high_modulated)); % Envolvente de la señal modulada

s_demodulada_low = filter(lowPassFilter, envelop_low);
s_demodulada_medium = filter(lowPassFilter, envelop_medium);
s_demodulada_high = filter(lowPassFilter, envelop_high);

s_normalizado_low = s_demodulada_low / max(abs(s_demodulada_low)); % Señal normalizada
s_normalizado_medium = s_demodulada_medium / max(abs(s_demodulada_medium)); % Señal normalizada
s_normalizado_high = s_demodulada_high / max(abs(s_demodulada_high)); % Señal normalizada

t_ruido_bajo = (0:length(s_normalizado_low)-1) / Fs;
t_ruido_medio = (0:length(s_normalizado_medium)-1) / Fs;
t_ruido_alto = (0:length(s_normalizado_high)-1) / Fs;

figure(3);
subplot(3,2,1)
plot(t_ruido_bajo,s_normalizado_low)
ylim([0.7 1]);
title('Señal demodulada con ruido bajo')

Nb = length(s_normalizado_low);                   
dftb = fftshift(fft(s_normalizado_low));          
f0b = (-Nb/2:Nb/2-1)*(Fs/Nb);      
DEPb = (1/(Fs*Nb))*abs(dftb).^2;  


% PSD Señal Demodulada con ruido bajo
subplot(3, 2, 2);
plot(f0b,10*log10(DEPb))
xlim([-1.5*fc 1.5*fc]);
ylim([-120 0]);
xlabel('f [Hz]');
ylabel('|S(f)|^2 [dB]');
title('PSD Señal Demodulada con ruido bajo')
grid on

subplot(3,2,3)
plot(t_ruido_medio,s_normalizado_medium)
ylim([0.7 1]);
title('Señal demodulada con ruido medio')


Nm = length(s_normalizado_medium);                   
dftm = fftshift(fft(s_normalizado_medium));          
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
plot(t_ruido_alto,s_normalizado_high)
ylim([0.6 1]);
title('Señal demodulada con ruido alto')

Nh = length(s_normalizado_high);                   
dfth = fftshift(fft(s_normalizado_high));          
f0h = (-Nh/2:Nh/2-1)*(Fs/Nh);      
DEPh = (1/(Fs*Nh))*abs(dfth).^2;  


% PSD Señal Demodulada con ruido alto
subplot(3, 2, 6);
plot(f0h,10*log10(DEPh))
xlim([-1.5*fc 1.5*fc]);
ylim([-120 0]);
xlabel('f [Hz]');
ylabel('|S(f)|^2 [dB]');
title('PSD Señal Demodulada con ruido alto')
grid on

figure (7)
subplot(4, 1, 1);
plot(t, s_normalizado );
ylim([0.7 1]);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Demodulada Sin Ruido');

subplot(4, 1, 2);
plot(t_ruido_bajo,s_normalizado_low );
ylim([0.7 1]);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Demodulada Con Ruido Bajo');

subplot(4, 1, 3);
plot(t_ruido_medio,s_normalizado_medium );
ylim([0.7 1]);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Demodulada Con Ruido Medio');

subplot(4, 1, 4);
plot(t_ruido_alto,s_normalizado_high );
ylim([0.6 1]);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Demodulada Con Ruido Alto');

%soundsc(s_normalizado,Fs)
%soundsc(s_normalizado_low,Fs)
%soundsc(s_normalizado_medium,Fs)
%soundsc(s_normalizado_high,Fs)