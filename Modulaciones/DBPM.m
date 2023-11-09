%%%%%%%%%%%%%%%%%%%%%%% WBPM (Wide Band Phase Modulation)%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

%Cargar la señal de audio

rootdirectory = 'Z:\Documents\TEC\II_semestre_2023\Comu I\Comu_1';
[m, fs] = audioread(fullfile(rootdirectory, 'AudioTarea.m4a'));

%Promediar los canales izquierdo y derecho
m = mean(m,2);
m = m(:);

% Parámetros de modulación
fc = 1000; % Frecuencia de la portadora en Hz
fs_mod = 44100; % Frecuencia de muestreo para la modulación
beta = 5; % Índice de modulación

% Modulación WBPM

Ac= 5;      % Amplitud de la señal portadora
t = (0:1/fs_mod:(length(m)-1)/fs_mod)';

c = Ac*cos(2*pi*fc*t');     % Señal portadora


%CON pmmod
y_mod = pmmod(m, fc, fs_mod, beta);

%Manualmente

%y_mod = c .* cos(2 * pi * fc * t + beta * m);


%y_demod = pmdemod(y_mod, fc, fs_mod, beta);
y_demod = pmdemod(y_mod, fc, fs_mod, beta);


subplot(3, 2, 3);
plot(t,c);
xlabel('Tiempo (s)');
xlim([6 6.025]);
ylabel('Amplitud');
ylim([-7 7]);
title('Señal portadora');

%Representacion de la portadora en la frecuencia
Nc = length(c);                   % Longitud de la señal c
dftc = fftshift(fft(c));          % Coloca la componente cero en el centro del espectro 
f02 = (-Nc/2:Nc/2-1)*(fs/Nc);      % Base de frecuencias centradas en 0
DEPc = (1/(fs*Nc))*abs(dftc).^2;  % Densidad espectral de potencia de c

subplot(3, 2, 4);
plot(f02,10*log10(DEPc))
xlim([-1.5*fc 1.5*fc]);
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
f01 = (-Nm/2:Nm/2-1)*(fs/Nm);      
DEPm = (1/(fs*Nm))*abs(dftm).^2;  

subplot(3, 2, 2);
plot(f01,10*log10(DEPm))
xlim([-1.5*fc 1.5*fc]);
ylim([-60 20]);
xlabel('f [Hz]');
ylabel('|M(f)|^2 [dB]');
title('Señal moduladora en la frecuencia ');
grid on

% Señal WBPM modulada en el tiempo
subplot(3, 2, 5);
plot(t, y_mod);
xlabel('Tiempo (s)');
xlim([6 6.025])
ylabel('Amplitud');
ylim([-2 2]);
title('Señal WBPM Modulada');

Ny = length(y_mod);                   
dfty = fftshift(fft(y_mod));          
f0 = (-Ny/2:Ny/2-1)*(fs/Ny);      
DEPy = (1/(fs*Ny))*abs(dfty).^2;  

% PSD de la señal WBPM modulada
subplot(3, 2, 6);
plot(f0,10*log10(DEPy))
xlim([-1.5*fc 1.5*fc]);
ylim([-50 10]);
xlabel('f [Hz]');
ylabel('|S(f)|^2 [dB]');
title('PSD señal DBPM modulada ')
grid on




%%%%%%%%%%%%%%%%%%%%%%%% Ruido %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Agregar ruido a las señales moduladas

snr_low = 30;       % Relación señal/ruido (bajo)
snr_medium = 15;     % Relación señal/ruido (medio)
snr_high = 8;       % Relación señal/ruido (alto)


% BWPM con ruido
y_mod_noisy_low = awgn(y_mod, snr_low);
y_mod_noisy_medium = awgn(y_mod, snr_medium);
y_mod_noisy_high = awgn(y_mod, snr_high);

%Modulación con ruido


%Demodulación con ruido 
y_demod_noisy_low = pmdemod(y_mod_noisy_low, fc, fs_mod, beta);
y_demod_noisy_medium = pmdemod(y_mod_noisy_medium, fc, fs_mod, beta);
y_demod_noisy_high = pmdemod(y_mod_noisy_high, fc, fs_mod, beta);

%Plot de señales demodulas

% Señal demodulada
figure(3);
subplot(4,1,1)
plot(t,y_demod)
title('Señal demodulada sin ruido')
xlabel('Tiempo (s)');
xlim([0 10])
ylabel('Amplitud');
ylim([-2 2]);

%Plot de demoduladoras con ruido
subplot(4,1,2)
plot(t,y_demod_noisy_low);
xlabel('Tiempo [s]');
xlim([0 10])
ylabel('Amplitud [V]');
ylim([-2 2])
title('Demoduladora con ruido bajo ')
subplot(4,1,3);
plot(t,y_demod_noisy_medium);
xlabel('Tiempo [s]');
xlim([0 10])
ylabel('Amplitud [V]');
ylim([-2 2])
title('Demoduladora con ruido medio ')
subplot(4,1,4);
plot(t,y_demod_noisy_high);
xlabel('Tiempo [s]');
xlim([0 10])
ylabel('Amplitud [V]');
ylim([-2 2])
title('Demoduladora con ruido alto ')



%Plots de modulación con ruido
figure(5);

subplot(4,1,1)
plot(t,y_mod);
xlabel('Tiempo [s]');
xlim([0 0.005])
ylabel('Amplitud [V]');
ylim([-2 2])
title('Modulación DBPM sin ruido')

subplot(4,1,2)
plot(t,y_mod_noisy_low);
xlabel('Tiempo [s]');
xlim([0 0.005])
ylabel('Amplitud [V]');
ylim([-2 2])
title('Modulación DBPM con ruido bajo')
subplot(4,1,3);
plot(t,y_mod_noisy_medium);
xlabel('Tiempo [s]');
xlim([0 0.005])
ylabel('Amplitud [V]');
ylim([-2 2])
title('Modulación DBPM con ruido medio')
subplot(4,1,4);
plot(t,y_mod_noisy_high);
xlabel('Tiempo [s]');
xlim([0 0.005])
ylabel('Amplitud [V]');
ylim([-2 2])
title('Modulación DBPM con ruido alto')



%Espectro sin ruido
Nc_no_noisy = length(y_demod);                        % Longitud de la señal c
dftc_no_noisy = fftshift(fft(y_demod));                         % Coloca la componente cero en el centro del espectro 
f02_no_noisy = (-Nc_no_noisy/2:Nc_no_noisy/2-1)*(fs/Nc_no_noisy);      % Base de frecuencias centradas en 0
DEPc_no_noisy = (1/(fs*Nc_no_noisy))*abs(dftc_no_noisy).^2;                       % Densidad espectral de potencia de c


%Espectro bajo ruido
Nc_noisy_low = length(y_demod_noisy_low);                        % Longitud de la señal c
dftc_noisy_low = fftshift(fft(y_demod_noisy_low));                         % Coloca la componente cero en el centro del espectro 
f02_noisy_low = (-Nc_noisy_low/2:Nc_noisy_low/2-1)*(fs/Nc_noisy_low);      % Base de frecuencias centradas en 0
DEPc_noisy_low = (1/(fs*Nc_noisy_low))*abs(dftc_noisy_low).^2;                       % Densidad espectral de potencia de c



%Espectro ruido medio
Nc_noisy_medium = length(y_demod_noisy_medium);                        % Longitud de la señal c
dftc_noisy_medium = fftshift(fft(y_demod_noisy_medium));                         % Coloca la componente cero en el centro del espectro 
f02_noisy_medium = (-Nc_noisy_medium/2:Nc_noisy_medium/2-1)*(fs/Nc_noisy_medium);      % Base de frecuencias centradas en 0
DEPc_noisy_medium = (1/(fs*Nc_noisy_medium))*abs(dftc_noisy_medium).^2;                       % Densidad espectral de potencia de c

%Espectro ruido alto

Nc_noisy_high = length(y_demod_noisy_high);                        % Longitud de la señal c
dftc_noisy_high = fftshift(fft(y_demod_noisy_high));                         % Coloca la componente cero en el centro del espectro 
f02_noisy_high = (-Nc_noisy_high/2:Nc_noisy_high/2-1)*(fs/Nc_noisy_high);      % Base de frecuencias centradas en 0
DEPc_noisy_high = (1/(fs*Nc_noisy_high))*abs(dftc_noisy_high).^2;   



% Espectros sin ruido
figure (8)
subplot(4,1,1)
plot(f02_no_noisy,10*log10(DEPc_no_noisy))
xlim([-1.5*fc 1.5*fc]);
ylim([-60 20]);
xlabel('f [Hz]');
ylabel('|R1(f)|^2 [dB]');
title('PSD señal demodulada sin ruido ')

%Espectros con ruidos
subplot(4,1,2)
plot(f02_noisy_low,10*log10(DEPc_noisy_low))
xlim([-1.5*fc 1.5*fc]);
ylim([-60 20]);
xlabel('f [Hz]');
ylabel('|R1(f)|^2 [dB]');
title('PSD señal demodulada con ruido bajo')
subplot(4,1,3)
plot(f02_noisy_medium,10*log10(DEPc_noisy_medium))
xlim([-1.5*fc 1.5*fc]);
ylim([-60 20]);
xlabel('f [Hz]');
ylabel('|R2(f)|^2 [dB]');
title('PSD señal demodulada con ruido medio')
subplot(4,1,4)
plot(f02_noisy_high,10*log10(DEPc_noisy_high))
xlim([-1.5*fc 1.5*fc]);
ylim([-60 20]);
xlabel('f [Hz]');
ylabel('|R3(f)|^2 [dB]');
title('PSD señal demodulada con ruido alto')


%Reproducir los audios demodulados con ruido

%sound(y_demod_noisy_low, fs);
%sound(y_demod_noisy_medium, fs);
%sound(y_demod_noisy_high, fs);




