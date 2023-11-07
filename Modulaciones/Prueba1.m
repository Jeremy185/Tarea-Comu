%%%%%%%%%%%%%%%%%%%%%%% NBFM (Banda Angosta) %%%%%%%%%%%%%%%%%%%%%%

%Cargar la señal de audio
clc;clear;
rootdirectory = 'Z:\Downloads';
[m, Fs] = audioread(fullfile(rootdirectory,'prueba2.m4a'));

m = m(:);

% Parámetros
Am = max(abs(m)); % Amplitud del mensaje
fm = 100; % Frecuencia del mensaje (Hz)
t = (0:length(m)-1) / Fs;


%Parametros de la señal portadora
fc = 1000;                     % Frecuencia de la señal portadora 1kHz
Ac = 1;                     % Amplitud de la señal portadora
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
f02 = (-Nc/2:Nc/2-1)*(Fs/Nc);      % Base de frecuencias centradas en 0
DEPc = (1/(Fs*Nc))*abs(dftc).^2;  % Densidad espectral de potencia de c


subplot(3, 2, 4);
plot(f02,10*log10(DEPc))
xlim([-1.5*fc 1.5*fc]);
ylim([-80 50]);
xlabel('f [Hz]');
ylabel('|C(f)|^2 [dB]');
title('Señal portadora en la frecuencia ');
grid on

% Representación de la moduladora en el tiempo

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
xlim([-1.5*fc 1.5*fc]);
ylim([-60 20]);
xlabel('f [Hz]');
ylabel('|M(f)|^2 [dB]');
title('Señal moduladora en la frecuencia ');
grid on

K=0.1;
%{
% Señal NBFM modulada en el tiempo
nbfm_modulated = c + K .*c1.*m;  %
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
title('PSD señal NBFM')
grid on

%Demodulacion%
y = pmdemod(nbfm_modulated, fc, Fs, K);
figure (2)
subplot(3, 1, 1);
plot(t, y);
xlabel('Tiempo (s)');
ylim([-1 1]);
ylabel('Amplitud');
title('Señal Demodulada');

sound (y, Fs)
%}


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

%}
% PSD de la señal NBFM modulada
subplot(3, 2, 6);
plot(f0,10*log10(DEPy))
xlim([-1.5*fc 1.5*fc]);
ylim([-120 0]);
xlabel('f [Hz]');
ylabel('|S(f)|^2 [dB]');
title('PSD señal NBFM')
grid on



%%%%%%%%%%%%%%%%%%%%%%%% Ruido %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Agregar ruido a las señales moduladas

snr_low = 5;       % Relación señal/ruido (bajo)
snr_medium = 1;     % Relación señal/ruido (medio)
snr_high = 0.5;       % Relación señal/ruido (alto)


% DSB-SC con ruido
nbfm_modulated_noisy_low = awgn(nbfm_modulated, snr_low);
nbfm_modulated_noisy_medium = awgn(nbfm_modulated, snr_medium);
nbfm_modulated_noisy_high = awgn(nbfm_modulated, snr_high);


figure(2);
subplot(3,1,1)
plot(t,nbfm_modulated_noisy_low);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación NBFM con ruido bajo ')
subplot(3,1,2);
plot(t,nbfm_modulated_noisy_medium);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación NBFM con ruido medio ')
subplot(3,1,3);
plot(t,nbfm_modulated_noisy_high);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación NBFM con ruido alto ')

%Demodulacion%
nbfm_demodulated_noisy_low = fmdemod(nbfm_modulated_noisy_low, fc, Fs, K);
nbfm_demodulated_noisy_medium = fmdemod(nbfm_modulated_noisy_medium, fc, Fs, K);
nbfm_demodulated_noisy_high = fmdemod(nbfm_modulated_noisy_high, fc, Fs, K);

% Demodulación con diferentes niveles de ruido %
figure(3);
subplot(3,2,1)
plot(t,nbfm_demodulated_noisy_low)
title('Señal demodulada con ruido bajo')
subplot(3,2,3)
plot(t,nbfm_demodulated_noisy_medium)
title('Señal demodulada con ruido medio')
subplot(3,2,5)
plot(t,nbfm_demodulated_noisy_high)
title('Señal demodulada con ruido alto')

%{
%%%%%%%%%%%%%%%%%%%%%%%% Ruido %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Agregar ruido a las señales moduladas

snr_low = 30;       % Relación señal/ruido (bajo)
snr_medium = 15;     % Relación señal/ruido (medio)
snr_high = 8;       % Relación señal/ruido (alto)


% DSB-SC con ruido
nbfm_modulated_noisy_low = awgn(nbfm_modulated, snr_low);
nbfm_modulated_noisy_medium = awgn(nbfm_modulated, snr_medium);
nbfm_modulated_noisy_high = awgn(nbfm_modulated, snr_high);

% Llamar a la función para aplicar el filtro pasa banda con Ruido bajo
[filtered_signal_r1, DEP_r1, f0_r1] = aplicarFiltroPasaBanda(nbfm_modulated_noisy_low, Fs, fc, fm);
% Llamar a la función para aplicar el filtro pasa banda con Ruido medio
[filtered_signal_r2, DEP_r2, f0_r2] = aplicarFiltroPasaBanda(nbfm_modulated_noisy_medium, Fs, fc, fm);
% Llamar a la función para aplicar el filtro pasa banda con Ruido alto
[filtered_signal_r3, DEP_r3, f0_r3] = aplicarFiltroPasaBanda(nbfm_modulated_noisy_high, Fs, fc, fm);

figure(2);
subplot(3,1,1)
plot(t,nbfm_modulated_noisy_low);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación NBFM con ruido bajo ')
subplot(3,1,2);
plot(t,nbfm_modulated_noisy_medium);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación NBFM con ruido medio ')
subplot(3,1,3);
plot(t,nbfm_modulated_noisy_high);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación NBFM con ruido alto ')

% Demodulación con diferentes niveles de ruido %
figure(3);
subplot(3,2,1)
plot(t,filtered_signal_r1)
title('Señal demodulada con ruido bajo')
subplot(3,2,3)
plot(t,filtered_signal_r2)
title('Señal demodulada con ruido medio')
subplot(3,2,5)
plot(t,filtered_signal_r3)
title('Señal demodulada con ruido alto')
% Espectros
subplot(3,2,2)
plot(f0_r1,10*log10(DEP_r1))
xlim([-1.5*fc 1.5*fc]);
ylim([-60 20]);
xlabel('f [Hz]');
ylabel('|R1(f)|^2 [dB]');
title('PSD señal demodulada con ruido bajo')
grid on

subplot(3,2,4)
plot(f0_r2,10*log10(DEP_r2))
xlim([-1.5*fc 1.5*fc]);
ylim([-60 20]);
xlabel('f [Hz]');
ylabel('|R2(f)|^2 [dB]');
title('PSD señal demodulada con ruido medio')
grid on

subplot(3,2,6)
plot(f0_r3,10*log10(DEP_r3))
xlim([-1.5*fc 1.5*fc]);
ylim([-60 20]);
xlabel('f [Hz]');
ylabel('|R3(f)|^2 [dB]');
title('PSD señal demodulada con ruido alto')
grid on

%Reproducir los audios demodulados con ruido

%soundsc(filtered_signal_r1, Fs);
%soundsc(filtered_signal_r2, Fs);
%soundsc(filtered_signal_r3, Fs);

%%%%%%%%%%%%%%%%%%%%%%%%% Demodulacion sin ruido aplicado %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%demodulated_signal = dsb_sc_modulated .* c;  % Demodulación usando la misma portadora
[filtered_signal, DEPz, f0z] = aplicarFiltroPasaBanda(nbfm_modulated, Fs, fc, fm);


% Señal demodulada
figure(4);
subplot(2,1,1)
plot(t,y)
title('Señal demodulada sin ruido')
grid on

% DEP de la señal demodulada
subplot(2,1,2)
plot(f0z,10*log10(DEPz))
xlim([-1.5*fc 1.5*fc]);
ylim([-50 20]);
xlabel('f [Hz]');ylabel('|Z(f)|^2 [dB]')
title('PSD señal demodulada sin ruido')
grid on

%Reproducir el audio demodulado
soundsc(filtered_signal, Fs);

%Demodulacion%
y = fmdemod(nbfm_modulated, fc, Fs, K);


%%%%%%%%%% Funcion para aplicar el filtro de demodulacion %%%%%%%%%%%%%%%%%

function [filtered_signal, DEP, f0] = aplicarFiltroPasaBanda(nbfm_modulated, Fs, fc, fm)
    % signal_to_demodulate: Señal a demodular
    % Fs: Frecuencia de muestreo
    % fc: Frecuencia central del filtro pasa banda
    % BW: Ancho de banda del filtro pasa banda
    K=0.5;
    signal_demodulated = fmdemod(nbfm_modulated, fc, Fs, K);
    BW = 2*fm;
    % Diseñar el filtro pasa banda
    h = fir1(90, [((fc - BW/2)/Fs), ((fc + BW/2)/Fs)]); %Diseña un filtro pasabanda de orden 90

    % Aplicar el filtro pasa banda a la señal
    filtered_signal = filter(h, 1, signal_demodulated);

    % Calcular la densidad espectral de potencia (DEP)
    Nz = length(filtered_signal);
    dftz = fftshift(fft(filtered_signal));
    f0 = (-Nz/2:Nz/2-1) * (Fs/Nz);
    DEP = (1/(Fs*Nz)) * abs(dftz).^2;
end

%}

