%%%%%%%%%%%%%%%%%%%%%%% DSB-SL (Portadora Suprimida) %%%%%%%%%%%%%%%%%%%%%%

%Cargar la señal de audio
%rootdirectory = 'C:\Users\jesus\OneDrive - Estudiantes ITCR\TEC\2023\II Semestre 2023\Comunicaciones Eléctricas I\Tarea\Modulacion';
[m, Fs] = audioread('Recording.mp3');

m = m(:);

% Parámetros

t = (0:length(m)-1) / Fs;
Am = max(abs(m)); % Amplitud máxima de la señal de audio.
fm = 216; %PARA ESTE AUDIO, frecuencia de la moduladora. 


%Parametros de la señal portadora
fc = 1000;                  % Frecuencia de la señal portadora 1kHz
Ac = 1;                     % Amplitud de la señal portadora
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
f02 = (-Nc/2:Nc/2-1)*(Fs/Nc);     % Base de frecuencias centradas en 0
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

% Señal DSB-SL modulada en el tiempo
Ka = 2;  %Sensitividad de modulacion
dsb_sl_modulated = (1+Ka*m).*c;       % s(t) = (1+ka*m(t))
sobre_mod = Ka*Am; 

kaux = 1/Am;% Valor maximo de Ka

if sobre_mod > 1
    disp('Hay sobremodulacion')
else 
    disp('Se esta modulando bien')
end 

subplot(3, 2, 5);
plot(t, dsb_sl_modulated);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal DSB-SL Modulada');

Ny = length(dsb_sl_modulated);                   
dfty = fftshift(fft(dsb_sl_modulated));          
f0 = (-Ny/2:Ny/2-1)*(Fs/Ny);      
DEPy = (1/(Fs*Ny))*abs(dfty).^2;  

% PSD de la señal DSB-SL modulada
subplot(3, 2, 6);
plot(f0,10*log10(DEPy))
xlim([-1.5*fc 1.5*fc]);
ylim([-50 20]);
xlabel('f [Hz]');
ylabel('|S(f)|^2 [dB]');
title('PSD señal DSB-SL modulada ')
grid on

%%%%%%%%%%%%%%%%%%%%%%%% Ruido %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Agregar ruido a las señales moduladas

snr_low = 30;       % Relación señal/ruido (bajo)
snr_medium = 15;     % Relación señal/ruido (medio)
snr_high = 8;       % Relación señal/ruido (alto)


% DSB-SL con ruido
dsb_sl_modulated_noisy_low = awgn(dsb_sl_modulated, snr_low);
dsb_sl_modulated_noisy_medium = awgn(dsb_sl_modulated, snr_medium);
dsb_sl_modulated_noisy_high = awgn(dsb_sl_modulated, snr_high);

% Llamar a la función para aplicar el filtro pasa banda con Ruido bajo
[filtered_signal_r1, DEP_r1, f0_r1] = aplicarFiltroPasaBanda(dsb_sl_modulated_noisy_low, Fs, fc, fm, c);
% Llamar a la función para aplicar el filtro pasa banda con Ruido medio
[filtered_signal_r2, DEP_r2, f0_r2] = aplicarFiltroPasaBanda(dsb_sl_modulated_noisy_medium, Fs, fc, fm, c);
% Llamar a la función para aplicar el filtro pasa banda con Ruido alto
[filtered_signal_r3, DEP_r3, f0_r3] = aplicarFiltroPasaBanda(dsb_sl_modulated_noisy_high, Fs, fc, fm, c);

figure(2);
subplot(3,1,1)
plot(t,dsb_sl_modulated_noisy_low);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación DSB-SL con ruido bajo ')
subplot(3,1,2);
plot(t,dsb_sl_modulated_noisy_medium);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación DSB-SL con ruido medio ')
subplot(3,1,3);
plot(t,dsb_sl_modulated_noisy_high);
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
title('Modulación DSB-SL con ruido alto ')

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

%soundsl(filtered_signal_r1, Fs);
%sound(filtered_signal_r2, Fs);
%soundsl(filtered_signal_r3, Fs);

%%%%%%%%%%%%%%%%%%%%%%%%% Demodulacion sin ruido aplicado %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[filtered_signal, DEPz, f0z] = aplicarFiltroPasaBanda(dsb_sl_modulated, Fs, fc, fm, c);


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
%sound(filtered_signal, Fs);


%%%%%%%%%% Funcion para aplicar el filtro de demodulacion %%%%%%%%%%%%%%%%%

function [filtered_signal, DEP, f0] = aplicarFiltroPasaBanda(signal_to_demodulate, Fs, fc, fm, c)
    % signal_to_demodulate: Señal a demodular
    % Fs: Frecuencia de muestreo
    % fc: Frecuencia central del filtro pasa banda
    % BW: Ancho de banda del filtro pasa banda
    % 90 es el orden del filtro, entre mas alto mas ideal

    signal_demodulated = signal_to_demodulate .* c;
    BW = 2*fm;

    % Diseñar el filtro pasa banda
    h = fir1(90, [((fc - BW/2)/Fs), ((fc + BW/2)/Fs)]);

    % Aplicar el filtro pasa banda a la señal
    filtered_signal = filter(h, 1, signal_demodulated);

    % Calcular la densidad espectral de potencia (DEP)
    Nz = length(filtered_signal);       
    dftz = fftshift(fft(filtered_signal)); %Transforma la señal a fourier (shift lo centra)
    f0 = (-Nz/2:Nz/2-1) * (Fs/Nz);         %Se crea un vector de frecuencia para la visualización.
    DEP = (1/(Fs*Nz)) * abs(dftz).^2;      %Elevar el valor absoluto al cuadrado
    %(1/(Fs*Nz)) es un factor de normalizacion para que la grafica sea mas
    %interpretable
end


