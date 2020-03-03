clc
clear all
close all

load('../resources/103m.mat');
load('../resources/labels.mat');
load('resampled_lpf.mat');
load('resampled_hpf.mat');
load('resampled_fir.mat');

Fs=360;
x_ECG(1,1:end) = x_ECG(1,1:end) - mean(x_ECG(1,1:end)); 

x_ECG360 = resample(x_ECG(1,:), 9, 5, resampled_fir, 1);

[z,p,k] = tf2zpk(resampled_fir,1);

figure('name','Resampled LPF ');
subplot(2,2,1);
zplane(z,p)
title('Poles and zeros diagram (Resampled LPF)');
grid on

x_n = zeros(1,Fs*3);
x_n(1) = 1;

y_n = filter(resampled_fir,1,x_n);
Y_w = fft(y_n);

subplot(2,2,2);
plot(linspace(-Fs/2,Fs/2,length(Y_w)),fftshift(abs(Y_w)))
title('Transfer function')
ylabel('Amplitude')
xlabel('Frecuency [Hz]')
grid on

subplot(2,2,3);
n=0:(length(y_n)-1);
stem(n,y_n(n+1),'MarkerFaceColor','blue')
title('Impulse response')
ylabel('Amplitude')
xlabel('n')
xlim([0 15])
grid on

subplot(2,2,4)
[gd_lp,w] = grpdelay(resampled_fir,1);
gd_lp = max(gd_lp)*ones(1,length(gd_lp));
plot(w/pi,gd_lp)
title('LPF Group delay')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Group delay (samples)')
grid on


[h,w]=freqz(resampled_fir,1);
[phi,w]=phasez(resampled_fir,1);

figure('name','Bode diagram (LPF)');

subplot(2,1,1);
plot(w/pi,20*log10(abs(h)))
title('Bode pasa-bajos (m�dulo)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude [dB]')
grid on

subplot(2,1,2);
plot(w/pi,phi)
title('Bode fase pasa-bajos')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Phase [rad]')
grid on

t = 0:1/Fs:(size(x_ECG360,2)-1)/Fs;
x_ECG = [x_ECG360 ; t];

marks = round(marks*9/5+gd_lp(end)/5);

%% Ejercicio 1
delay_factor=0;
samples=12000;
m_inf = delay_factor+1; 
m_sup = delay_factor+round(samples*9/5);
xlim_sup = m_sup/Fs;
xlim_inf = m_inf/Fs;

figure('name', 'ECG manually labeled')
plot(x_ECG(2,m_inf:m_sup),x_ECG(1,m_inf:m_sup))
title('ECC manually labeled');
xlabel('Time [s]');
ylabel('Voltage [mV]');
grid on
hold on
plot(x_ECG(2,round(marks)),x_ECG(1,round(marks)),'r x','MarkerSize',6)
xlim([xlim_inf xlim_sup]);  
 

b1 = 36.0164 * resampled_lpf;
a1 = 1
x = filter(b1,a1,x_ECG(1,m_inf:m_sup));

figure('name','Se�al filtrada LP');
plot(x_ECG(2,m_inf:m_sup),x)
title('Se�al filtrada con el lowpass filter');
xlabel('Time [s]');
ylabel('Voltage [mV]');
xlim([xlim_inf xlim_sup]);
grid on

[z,p,k] = tf2zpk(b1,a1);

figure('name','Caracter�sticas del pasa-bajos');

subplot(2,2,1);
zplane(z,p)
title('Poles and zeros diagram (lowpass filter)');
grid on

x_n = zeros(1,Fs*3);
x_n(1) = 1;

y_n = filter(b1,a1,x_n);
Y_w = fft(y_n);

subplot(2,2,2);
plot(linspace(-Fs/2,Fs/2,length(Y_w)),fftshift(abs(Y_w)))
title('Transfer function')
ylabel('Amplitude')
xlabel('Frecuency [Hz]')
xlim([-Fs/2 Fs/2])
grid on

subplot(2,2,3);
n=0:(length(y_n)-1);
stem(n,y_n(n+1),'MarkerFaceColor','blue')
title('Impulse response')
ylabel('Amplitude')
xlabel('n')
xlim([0 15])
grid on

subplot(2,2,4)

[gd_lp,w] = grpdelay(b1,a1);
gd_lp = max(gd_lp)*ones(1,length(gd_lp)); %problemas con ploteo
plot(w/pi,gd_lp)
title('Group delay factor pasa-bajos')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Group delay (samples)')
grid on

% By looking at the Bode diagram it's possible to spot the -3 dB
% point and Fc is calculated performing f*Fs/2.
% fc = 19.33 Hz

figure('name','LPF Bode Diagram');

[h, w] = freqz(b1,a1);

subplot(2,1,1);
plot(w/pi,20*log10(abs(h)))
title('LPF Bode diagram (magnitude)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude [dB]')
grid on

[phi, w] = phasez(b1,a1);

subplot(2,1,2);
plot(w/pi,phi)
title('LPF Bode diagram (phase)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Phase [rad]')
grid on

b2 = resampled_hpf;
a2 = 1;

y = filter(b2,a2,x);

figure('name','HPF Output');

plot(x_ECG(2,m_inf:m_sup),y)
title('HPF Output');
xlabel('Time [s]');
ylabel('Voltage [mV]');
xlim([xlim_inf xlim_sup]);
grid on

[z,p,k] = tf2zpk(b2,a2); 

z_aux=zeros(31,1);
p_aux=zeros(31,1);

z_aux(1:16) = z(1:16);
z_aux(17:31)= z(18:32);
p_aux = p(1:31);

figure('name','HPF Characteristics');
subplot(2,2,1);
zplane(z_aux,p_aux)
title('Poles and zeros diagram (HPF)');
grid on

x_n = zeros(1,Fs*3);
x_n(1) = 1;

y_n = filter(b2,a2,x_n);
Y_w = fft(y_n);

subplot(2,2,2);
plot(linspace(-Fs/2,Fs/2,length(Y_w)),fftshift(abs(Y_w)))
title('Transfer function')
ylabel('Amplitude')
xlabel('Frecuency [Hz]')
xlim([-Fs/2 Fs/2])
grid on

subplot(2,2,3)
n=0:(length(y_n)-1);
stem(n,y_n(n+1), 'MarkerFaceColor','blue')
title('Impulse response')
ylabel('Amplitude')
xlabel('n')
xlim([0 40])
grid on

subplot(2,2,4)
[gd_hp,w] = grpdelay(b2,a2);

plot(w/pi,gd_hp)
title('Group delay factor pasa-altos')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Group delay (samples)')
grid on

% By looking at the Bode diagram it's possible to spot the -3 dB
% point and Fc is calculated performing f*Fs/2.
% fc = 8.42 Hz

figure('name','HPF Bode Diagram')
subplot(2,1,1);

[h, w] = freqz(b2,a2);

plot(w/pi,20*log10(abs(h)))
title('HPF Bode diagram (magnitude)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude [dB]')
grid on

[phi, w] = phasez(b2,a2);

subplot(2,1,2);
plot(w/pi,phi)
title('HPF Bode diagram (phase)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Phase [rad]')
grid on

bd = (Fs/8) * [1 2 0 -2 -1];
ad = 1;

x_dot = filter(bd,ad,y);
figure('name','Differentiator Output');
plot(x_ECG(2,m_inf:m_sup),x_dot)
title('Differentiator Output');
xlabel('Time [s]');
ylabel('Amplitude');
xlim([xlim_inf xlim_sup]);
grid on

[z,p,k] = tf2zpk(bd,ad);
x_n = zeros(1,Fs*3);
x_n(1) = 1;

y_n = filter(bd,ad,x_n);
Y_w = fft(y_n);

figure('name','Differentiator Characteristics')
subplot(2,2,1);
zplane(z,p)
title('Poles and zeros diagram (differentiator)');

subplot(2,2,2);
plot(linspace(-Fs/2,Fs/2,length(Y_w)),fftshift(abs(Y_w)))
title('Transfer function')
ylabel('Amplitude')
xlabel('Frecuency [Hz]')

grid on

subplot(2,2,3)
n=0:(length(y_n)-1);
stem(n,y_n(n+1), 'MarkerFaceColor','blue')
title('Impulse response')
ylabel('Amplitude')
xlabel('n')
xlim([0 40])
grid on

subplot(2,2,4)
[gd_d,w] = grpdelay(bd,ad);

plot(w/pi,gd_d)
title('Group delay differentiator')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Group delay (samples)')
grid on

figure('name','Diagrama de Bode del differentiator');

[h,w]=freqz(bd,ad);

subplot(2,1,1)
plot(w/pi,20*log10(abs(h)))
title('M�dulo differentiator')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude [dB]')
grid on

[phi,w]=phasez(bd,ad);

subplot(2,1,2)
plot(w/pi,phi)
title('Phase [rad]')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Phase [rad]')
grid on

x_sqr = (x_dot).^2;
X_sqr = fft(x_sqr); 

figure('name','Se�al al cuadrado');
subplot(2,1,1)
plot(x_ECG(2,m_inf:m_sup),x_sqr);
title('Squared Differentiator Output');
xlabel('Time [s]');
ylabel('Voltage [mV]');
xlim([xlim_inf xlim_sup]);
grid on 

subplot(2,1,2)
stem(linspace(-Fs/2,Fs/2,length(X_sqr)),fftshift(abs(X_sqr)),'MarkerFaceColor','blue')
title('Squared signal Spectrum');
ylabel('Amplitude');
xlabel('Frecuency [Hz]');
grid on

figure('name','Integrator Output');

n_step = round(31*9/5);
for n_integ = n_step*4:-n_step:n_step
  bi = 1/n_integ*ones(1,n_integ);
  ai = 1;
  x_integ = filter(bi,ai,x_sqr);
  subplot(2,2,n_integ/n_step)
  plot(x_ECG(2,m_inf:m_sup),x_integ)
  title( ['Integrator output N = ' num2str( n_integ )])
  xlabel('Time [s]');
  ylabel('Amplitude');
  xlim([xlim_inf xlim_sup]);
  grid on  
end

[z,p,k] = tf2zpk(bi,ai);
x_n = zeros(1,Fs*3);
x_n(1) = 1;

y_n = filter(bi,ai,x_n);
Y_w = fft(y_n);

figure('name','Integrator Characteristics')
subplot(2,2,1);
zplane(z,p)
title('Poles and zeros diagram (differentiator)');

subplot(2,2,2);
plot(linspace(-Fs/2,Fs/2,length(Y_w)),fftshift(abs(Y_w)))
title('Transfer function')
ylabel('Amplitude')
xlabel('Frecuency [Hz]')
grid on

subplot(2,2,3)
n=0:(length(y_n)-1);
stem(n,y_n(n+1), 'MarkerFaceColor','blue')
title('Impulse response')
ylabel('Amplitude')
xlabel('n')
xlim([0 40])
grid on

subplot(2,2,4)
[gd_i,w] = grpdelay(bi,ai);

plot(w/pi,gd_i)
title('Integrator Group delay')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Group delay (samples)')
grid on

[h, w] = freqz(bi,ai);
[phi, w] = phasez(bi,ai);

figure('name','Integrator Bode diagram');

subplot(2,1,1)
plot(w/pi,20*log10(abs(h)))
title('Integrator Bode diagram (module)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude [dB]')  
grid on

subplot(2,1,2)
plot(w/pi,phi)
title('Integrator Bode diagram (phase)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Phase (radians)')  
grid on
