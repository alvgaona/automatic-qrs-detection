clc
clear
close all

load('../resources/103m.mat');
load('../resources/labels.mat');

% Previously the 103m.mat file had only ECG values in millivolts.
% The signal was sampled at 200 Hz and another row was added with
% equally spaced samples 0.005 ms.

% The file 103m.mat is a 2x361112 matrix in which the first row is the
% amplitude of the signal and the second row the times each sample was
% recorded.

Fs = 200;             % Sampling frequency
num_samples = 361112; % Number of samples

x(1,1:end) = x(1,1:end) - mean(x(1,1:end)); 

%% Rubric point: Manual QRS complex detection

delay_factor = 12000;
samples = 12000;
m_inf = delay_factor + 1; 
m_sup = delay_factor + samples;

% Plot limits in time axis
xlim_sup = m_sup / Fs;
xlim_inf = m_inf / Fs;

fftpoints = m_sup + 1 - m_inf;

figure('name','ECG Signal');
plot(x(2,m_inf:m_sup),x(1,m_inf:m_sup))
title('ECG signal with manually labeled QRS complexes');
xlabel('Time [s]');
xlim([xlim_inf xlim_sup]);
ylabel('Voltage [mV]');
grid on
hold on

plot(x(2,ceil(labels)),x(1,ceil(labels)),'r x','MarkerSize',6)

%% Rubric point: QRS Frequency characteristics

qrs = x(1,5705:5720); % One QRS complex
QRS = fft(qrs, fftpoints);

figure('name','FFT abs(QRS)');
plot(linspace(-Fs/2,Fs/2,length(QRS)),fftshift(abs(QRS)))
title('QRS spectrum');
ylabel('Amplitude');
xlabel('Frequency [Hz]');
grid on

figure('name','FFT phase(QRS)');
plot(linspace(-Fs/2,Fs/2,length(QRS)),fftshift(angle(QRS)))
title('QRS spectrum');
ylabel('Phase');
xlabel('Frequency [Hz]');
grid on

%% Rubric point: ECG Spectrogram

figure('name','Spectrogram');

% Plotting whole signal including noise between 1200 and 1600 seconds.
subplot(2,2,1)
N = 2^6;
win = window(@hamming, N); 
spectrogram(x(1,:), win, N/2, size(x,1), Fs, 'yaxis');
xlabel('Time [s]');
ylabel('Frequency (Hz)');
title('Signal Spectrogram N=2^6')
grid on

subplot(2,2,2)
N = 2^7;
win = window(@hamming, N);                                                             
spectrogram(x(1,m_inf:m_sup), win, N/2, fftpoints, Fs, 'yaxis');
xlabel('Time [s]');
ylabel('Frequency (Hz)');
title('Signal Spectrogram N=2^7')
grid on

subplot(2,2,3)
N = 2^6;
win = window(@hamming, N); 
spectrogram(x(1,m_inf:m_sup), win, N/2, fftpoints, Fs, 'yaxis');
xlabel('Time [s]');
ylabel('Frequency (Hz)');
title('Signal Spectrogram N=2^6')
grid on

subplot(2,2,4)
N = 2^10;
win = window(@hamming, N); 
spectrogram(x(1,m_inf:m_sup), win, N/2, fftpoints, Fs, 'yaxis');
xlabel('Time [s]');
ylabel('Frequency (Hz)');
title('Signal Spectrogram N=2^10')
grid on

%% Rubric point: ECG Noise

figure('name','ECG Spectrum');
X = fft(x(1,m_inf:m_sup));
stem(linspace(-Fs/2,Fs/2,length(X)),fftshift(abs(X)),'-ko','MarkerSize',2);
title('ECG Spectrum');
ylabel('Amplitude');
xlabel('Frequency [Hz]');
grid on

figure('name','ECG Signal (temporal representation)');
plot(x(2,:), x(1,:))
title('ECG Signal');
ylabel('Amplitude');
xlabel('Time [s]');
grid on

%% Rubric point: Low-pass Filter

b1 = [1 0 0 0 0 0 -2 0 0 0 0 0 1]; % LPF numerator
a1 = [1 -2 1];                     % LPF denominator

y1 = filter(b1,a1,x(1,:)); % Filtrado

figure('name','Se�al filtrada LP');
plot(x(2,m_inf:m_sup),y1(m_inf:m_sup))
title('Low-pass filter output');
xlabel('Time [s]');
ylabel('Voltage [mV]');
xlim([xlim_inf xlim_sup]);
grid on

[z,p,~] = tf2zpk(b1,a1);

z_aux = zeros(10,1);
p_aux = zeros(10,1);

z_aux(1:6) = z(1:6);
z_aux(7:10)= z(9:12);
p_aux = p(1:10);

figure('name','LPF characteristics');

subplot(2,2,1);
zplane(z_aux,p_aux)
title('Pole and zeros diagram (LPF)');
grid on

x_n = zeros(1,Fs*3);
x_n(1) = 1;

y_n = filter(b1,a1,x_n);
Y_w = fft(y_n);

subplot(2,2,2);
plot(linspace(-Fs/2,Fs/2,length(Y_w)),fftshift(abs(Y_w)))
title('LPF Spectrum H(f)')
ylabel('Amplitude')
xlabel('Frequency [Hz]')
grid on

subplot(2,2,3);
n = 0:(length(y_n) - 1);
stem(n,y_n(n+1),'MarkerFaceColor','blue')
title('Impulse response h[n]')
ylabel('Amplitud')
xlabel('n')
xlim([0 15])
grid on

subplot(2,2,4)
[gd_lp, w] = grpdelay(b1,a1);
gd_lp = max(gd_lp)*ones(1,length(gd_lp));
plot(w/pi,gd_lp)
title('LPF group delay')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Group delay')
grid on

% By looking at the Bode diagram it's possible to spot the -3 dB
% point and Fc is calculated performing f*Fs/2.
% Fc = 10.74 Hz

figure('name','LPF Bode diagram');

[h,~] = freqz(b1,a1);
[phi,w] = phasez(b1,a1);

subplot(2,1,1);
plot(w/pi,20*log10(abs(h)))
title('LPF Bode diagram (magnitude)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
grid on

subplot(2,1,2);
plot(w/pi,phi)
title('LPF Bode diagram (phase)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Phase [rad]')
grid on

%% Rubric point: High-pass filter

b2 = [-1/32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1/32]; % HPF Numerator
a2 = [1 -1];                                                                      % HPF Denomintaor

y2 = filter(b2,a2,y1); 

figure('name','HPF Output');

plot(x(2,m_inf:m_sup),y2(m_inf:m_sup))
title('High-pass filter output');
xlabel('Time [s]');
ylabel('Voltage [mV]');
xlim([xlim_inf xlim_sup]);
grid on

[z,p,~] = tf2zpk(b2,a2); 

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
title('HPF Spectrum - H(f)')
ylabel('Amplitude')
xlabel('Frequency [Hz]')
grid on

subplot(2,2,3)
n=0:(length(y_n)-1);
stem(n,y_n(n+1), 'MarkerFaceColor','blue')
title('Impulse response h[n])')
ylabel('Amplitude')
xlabel('Time (n)')
xlim([0 40])
grid on

subplot(2,2,4)
[gd_hp,w] = grpdelay(b2,a2);

plot(w/pi,gd_hp)
title('HPF Group Delay')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Group delay (samples)')
grid on

% By looking at the Bode diagram it's possible to spot the -3 dB
% point and Fc is calculated performing f*Fs/2.
% Fc = 4.68 Hz

figure('name','Diagrama de Bode del pasa-altos')
subplot(2,1,1);

[h,w]=freqz(b2,a2);

plot(w/pi,20*log10(abs(h)))
title('HPF Spectrum (magnitude)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
grid on

[phi,w]=phasez(b2,a2);

subplot(2,1,2);
plot(w/pi,phi)
title('HPF Spectrum (phase)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Phase [rad]')
grid on

%% Rubric point: FIR Filters

% Polynomes division performed with devconv
LPFCoeff = deconv(b1,a1); 
HPFCoeff = deconv(b2,a2);

x_FIRLP = filter(LPFCoeff,1,x(1,m_inf:m_sup));
x_FIRHP = filter(HPFCoeff,1,x_FIRLP);

figure('name','FIR LPF Output');
plot(x(2,m_inf:m_sup),x_FIRLP)
title('FIR LPF Output');
xlabel('Time [s]');
ylabel('Voltage [mV]');
xlim([xlim_inf xlim_sup]);
grid on

figure('name','FIR HPF Output');
plot(x(2,m_inf:m_sup),x_FIRHP)
title('FIR HPF Output');
xlabel('Time [s]');
ylabel('Voltage [mV]');
xlim([xlim_inf xlim_sup]);
grid on

%% Rubric point: ECG derivative

bd = (Fs/8) * [1 2 0 -2 -1];
ad = 1;

x_dot = filter(bd, ad, y2);
figure('name','Differentiator Output');
plot(x(2,m_inf:m_sup), x_dot(m_inf:m_sup))
title('Differentiator Output');
xlabel('Time [s]');
ylabel('Amplitude');
xlim([xlim_inf xlim_sup]);
grid on

[z,p,~] = tf2zpk(bd,ad);
x_n = zeros(1,Fs*3);
x_n(1) = 1;

y_n = filter(bd,ad,x_n);
Y_w = fft(y_n);

figure('name','Differentatior Characteristics')
subplot(2,2,1);
zplane(z,p)
title('Poles and zeros diagram (differentiator)');

subplot(2,2,2);
plot(linspace(-Fs/2,Fs/2,length(Y_w)),fftshift(abs(Y_w)))
title('Transfer function Spectrum - H(f)')
ylabel('Amplitude')
xlabel('Frequency [Hz]')
grid on

subplot(2,2,3)
n = 0:(length(y_n)-1);
stem(n, y_n(n+1), 'MarkerFaceColor','blue')
title('Impulse response - h[n])')
ylabel('Amplitude')
xlabel('n')
xlim([0 40])
grid on

subplot(2,2,4)
[gd_d,w] = grpdelay(bd,ad);
plot(w/pi,gd_d)
title('Differentiator Group Delay')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Group delay (samples)')
grid on

[h,~]=freqz(bd,ad);
[phi,w]=phasez(bd,ad);

figure('name','Differentiator Bode Diagram');
subplot(2,1,1)
plot(w/pi,20*log10(abs(h)))
title('Differentiator Bode Diagram (magnitude)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
grid on

subplot(2,1,2)
plot(w/pi,phi)
title('Differentiator Bode Diagram (phase)')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Phase [rad]')
grid on

%% Rubric point: Squared Output
x_sqr = (x_dot).^2;
X_sqr = fft(x_sqr); 

figure('name','Squared Output');
subplot(2,1,1)
plot(x(2,m_inf:m_sup),x_sqr(m_inf:m_sup));
title('Squared Output');
xlabel('Time [s]');
ylabel('Voltage [mV]');
xlim([xlim_inf xlim_sup]);
grid on 

subplot(2,1,2)
stem(linspace(-Fs/2,Fs/2,length(X_sqr)),fftshift(abs(X_sqr)),'MarkerFaceColor','blue')
title('Squared Signal Spectrum');
ylabel('Amplitude');
xlabel('Frequency [Hz]');
grid on

%% Rubric point: Integrator System Frequency Response

figure('name','Integrator Output');

n_step = 31;

for n_integ = n_step*4:-n_step:n_step
  bi = 1/n_integ*ones(1,n_integ);
  ai = 1;
  x_integ = filter(bi,ai,x_sqr); %promedio ventaneado
  subplot(2,2,n_integ/n_step)
  plot(x(2,m_inf:m_sup),x_integ(m_inf:m_sup))
  title( ['Signal through the integrator N=' num2str( n_integ )])
  xlabel('Time [s]');
  ylabel('Amplitude');
  xlim([xlim_inf xlim_sup]);
  grid on  
end

% Poles and zeros
[z,p,k] = tf2zpk(bi,ai);
x_n = zeros(1,Fs*3);
x_n(1) = 1;

y_n = filter(bi,ai,x_n);
Y_w = fft(y_n);

figure('name','Integrator Characteristics')
subplot(2,2,1);
zplane(z,p)
title('Poles and zeros diagrama (integrator)');

subplot(2,2,2);
plot(linspace(-Fs/2,Fs/2,length(Y_w)),fftshift(abs(Y_w)))
title('Transfer function Spectrum - H(f)')
ylabel('Amplitude')
xlabel('Frequency (f)[Hz]')
grid on

subplot(2,2,3)
n = 0:(length(y_n)-1);
stem(n,y_n(n+1), 'MarkerFaceColor','blue')
title('Impulse response - h[n])')
ylabel('Amplitude')
xlabel('n')
xlim([0 40])
grid on

subplot(2,2,4)
[gd_i,w] = grpdelay(bi,ai);

plot(w/pi,gd_i)
title('Integrator group delay')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Group delay (samples)')
grid on


[h,~]=freqz(bi,ai);
[phi,w] = phasez(bi,ai);

figure('name','Bode integrador');

subplot(2,1,1)
plot(w/pi,20*log10(abs(h)))
title('Bode m�dulo integrador')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')  
grid on

subplot(2,1,2)
plot(w/pi,phi)
title('Bode fase integrador')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Phase (radians)')  
grid on

%% Rubric point: Automatic Detection

labels = rwavedet(y2,2000,5);

figure('name','R Wave Detector')
plot(x(2,m_inf:m_sup),y2(m_inf:m_sup))
title('Signal with labeled peaks');
xlabel('Time [s]');
ylabel('Voltage [mV]');
grid on
xlim([xlim_inf xlim_sup]);
hold on

labels_pos = find(labels);

plot(x(2,labels_pos),y2(labels_pos),'mo')

%% Rubric point: Gold Standard

bandpass_grpd = gd_lp(end) + gd_hp(end);
bandpass_grpd = ceil(bandpass_grpd);

aux = (labels>m_inf & labels<m_sup);
new_labels = find(aux);  

estimated_marks = pan_tompkin(x_integ, y2, Fs);

figure('name','Gold Standard Automatic Detection');
hold on
plot(x(2,m_inf:m_sup), y2(m_inf:m_sup));
xlim([xlim_inf xlim_sup]);
plot(x(2,estimated_marks), y2(estimated_marks), 'r x');
plot(x(2,round(labels(new_labels)+bandpass_grpd)), y2(round(labels(new_labels) + bandpass_grpd)),'g o', 'MarkerSize', 6);
title('Signal with peaks labeled');
legend('Location','Southeast','Filtered ECG','Gold Standard','Manual labels')
xlabel('Time [s]');
ylabel('Voltage [mV]');
grid on

% Performance comparison
[TP, FP, FN] = performance(round(labels) + bandpass_grpd, labels_pos, 5);
[TP_GS, FP_GS, FN_GS] = performance(round(labels) + bandpass_grpd, estimated_marks, 5);

%% Rubric point: Signal-to-noise Ratio

% Noise added: SNR 10 dB, 20 dB, 30 dB

SNRi = 10; % Initial SNR
figure('name','Gold Standard detection with added noise.');
x_reference=y2;

for i = 1:3
  noise = randn(1,length(x(2,:)));
  SNRd = SNRi * i; % Desired SNR
  
  SNR=10*log10(sum(x(1,:).^2)/sum(noise(1,:).^2));
  cte=sqrt(10^((SNR-SNRd ) /10));

  noise=cte*noise;
  x_noise=x(1,:)+noise(1,:);
  
  y = filter(b1,a1,x_noise);
  y2 = filter(b2,a2,y);
  x_dot = filter(bd,ad,y2);
  x_sqr = (x_dot).^2;
  x_integ = filter(bi,ai,x_sqr);

  noise_filter1 = filter(b1,a1,noise(1,:));
  noise_filter2 = filter(b2,a2,noise_filter1);

  SNRf=10*log10(sum(x_reference.^2)/sum(noise_filter2.^2));
  estimated_marks = pan_tompkin(x_integ,y2,Fs);
  [TP_noise(i), FP_noise(i), FN_noise(i)] = performance(round(labels) + bandpass_grpd, estimated_marks, 5);

  subplot(3,1,i)  
  plot(x(2,m_inf:m_sup), y2(m_inf:m_sup));
  hold on
  plot(x(2,estimated_marks), y2(estimated_marks), 'r x');
  plot(x(2,round(labels(new_labels) + bandpass_grpd)), y2(round(labels(new_labels) + bandpass_grpd)),'g o', 'MarkerSize', 6);
  str = sprintf('SNR_{desired} = %d dB', SNRd) ;
  title(str);
  xlabel('Time [s]');
  ylabel('Voltage [mV]');
  xlim([xlim_inf xlim_sup]);
  grid on
  legend('Location','Southeast','Se�al ECG filtrada','Gold Standard','labels manuales')
  str = sprintf('SNR_{final} = %d dB', round(SNRf));
  text(xlim_sup,6000,str)
  str = sprintf('FP = %d', fp_noise(i));
  text(xlim_sup,4000,str)
  str = sprintf('FN = %d', fn_noise(i));
  text(xlim_sup,2000,str)
end
