%%%%%%%%%%%%%%%%%%%%%%% DSB-SC (Portadora Suprimida) %%%%%%%%%%%%%%%%%%%%%%

%Cargar la señal de audio
rootdirectory = 'C:\Users\jesus\OneDrive - Estudiantes ITCR\TEC\2023\II Semestre 2023\Comunicaciones Eléctricas I\Tarea\Modulacion';
[m, Fs] = audioread(fullfile(rootdirectory, 'prueba.m4a'));

m = m(:);

% Parámetros
t = (0:length(m)-1) / Fs;
Am = max(abs(m)); % Amplitud máxima de la señal de audio
fm = 120; %PARA ESTE AUDIO

%Parametros de la señal portadora
fc = 1000;                  % Frecuencia de la señal portadora 1kHz
Ac = 8;                     % Amplitud de la señal portadora
c = Ac*cos(2*pi*fc*t');     % Señal portadora
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
xlim([-1.5*fc 1.5*fc]);
ylim([-60 20]);
xlabel('f [Hz]');
ylabel('|M(f)|^2 [dB]');
title('Señal moduladora en la frecuencia ');
grid on

% Señal DSB-SC modulada en el tiempo
dsb_sc_modulated = m .*c;       % s(t) = m(t) * c(t)
subplot(3, 2, 5);
plot(t, dsb_sc_modulated);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal DSB-SC Modulada');

Ny = length(dsb_sc_modulated);                   
dfty = fftshift(fft(dsb_sc_modulated));          
f0 = (-Ny/2:Ny/2-1)*(Fs/Ny);      
DEPy = (1/(Fs*Ny))*abs(dfty).^2;  

% PSD de la señal DSB-SC modulada
subplot(3, 2, 6);
plot(f0,10*log10(DEPy))
xlim([-1.5*fc 1.5*fc]);
ylim([-50 20]);
xlabel('f [Hz]');
ylabel('|S(f)|^2 [dB]');
title('PSD señal DSB-SC modulada ')
grid on

%%%%%%%%%%%%%%%%%%%%%%%% Ruido %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Agregar ruido a las señales moduladas

snr_low = 30;       % Relación señal/ruido (bajo)
snr_medium = 15;     % Relación señal/ruido (medio)
snr_high = 8;       % Relación señal/ruido (alto)


% DSB-SC con ruido
dsb_sc_modulated_noisy_low = awgn(dsb_sc_modulated, snr_low);
dsb_sc_modulated_noisy_medium = awgn(dsb_sc_modulated, snr_medium);
dsb_sc_modulated_noisy_high = awgn(dsb_sc_modulated, snr_high);

% Llamar a la función para aplicar el filtro pasa banda con Ruido bajo
[filtered_signal_r1, DEP_r1, f0_r1] = aplicarFiltroPasaBanda(dsb_sc_modulated_noisy_low, Fs, fc, fm,c);
% Llamar a la función para aplicar el filtro pasa banda con Ruido medio
[filtered_signal_r2, DEP_r2, f0_r2] = aplicarFiltroPasaBanda(dsb_sc_modulated_noisy_medium, Fs, fc, fm,c);
% Llamar a la función para aplicar el filtro pasa banda con Ruido alto
[filtered_signal_r3, DEP_r3, f0_r3] = aplicarFiltroPasaBanda(dsb_sc_modulated_noisy_high, Fs, fc, fm,c);

figure(2);
subplot(3,1,1)
plot(t,dsb_sc_modulated_noisy_low);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación DSB-SC con ruido bajo ')
subplot(3,1,2);
plot(t,dsb_sc_modulated_noisy_medium);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación DSB-SC con ruido medio ')
subplot(3,1,3);
plot(t,dsb_sc_modulated_noisy_high);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación DSB-SC con ruido alto ')

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
[filtered_signal, DEPz, f0z] = aplicarFiltroPasaBanda(dsb_sc_modulated, Fs, fc, fm,c);


% Señal demodulada
figure(4);
subplot(2,1,1)
plot(t,filtered_signal)
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
%soundsc(filtered_signal, Fs);


%%%%%%%%%% Funcion para aplicar el filtro de demodulacion %%%%%%%%%%%%%%%%%

function [filtered_signal, DEP, f0] = aplicarFiltroPasaBanda(signal_to_demodulate, Fs, fc, fm,c)
    % signal_to_demodulate: Señal a demodular
    % Fs: Frecuencia de muestreo
    % fc: Frecuencia central del filtro pasa banda
    % BW: Ancho de banda del filtro pasa banda
    sginal_modulated = signal_to_demodulate .* c;
    BW = 2*fm;
    % Diseñar el filtro pasa banda
    h = fir1(90, [((fc - BW/2)/Fs), ((fc + BW/2)/Fs)]); %Diseña un filtro pasabanda de orden 90

    % Aplicar el filtro pasa banda a la señal
    filtered_signal = filter(h, 1, sginal_modulated);

    % Calcular la densidad espectral de potencia (DEP)
    Nz = length(filtered_signal);
    dftz = fftshift(fft(filtered_signal));
    f0 = (-Nz/2:Nz/2-1) * (Fs/Nz);
    DEP = (1/(Fs*Nz)) * abs(dftz).^2;
end