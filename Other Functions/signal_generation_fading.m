% This script deals as signal generation for Fading channels

clear all;
close all;
tic
%% Define main parameters

% Carrier Frequency
fc = 55e3;
% Sampling Frequency
fs = 150e3;
% Number of Samples
S = 5e7;
% Number of samples of the window with constant fading factor
N = 5e3;
% Time vector (discrete)
n = (0:1:S-1);
% Noise Floor
sigma_w = sqrt(7.5e-10);
% Scale factor of Fading
h_sigma = 1; % Can it be lower than 1? Whats the meaning?

% Obtain amplitudes depending on the SNR (-5, 0, ..., 20)
A = zeros(1,6);
for i = 0:5
    A(i+1) = sqrt(sigma_w^2*10^(0.1*(-5+5*i)));
end

tic
x_m5 = 0.05*exp(2i*pi*fc/fs*n);
h_m5 = h_sigma*(randn(1,S/N)+1i*randn(1,S/N));
h_m5 = repelem(h_m5,N);
s_m5 = h_m5.*x_m5; 
% save('Files/Rayleigh Fading/1 Canal/Enviat/final_rayfad_1','s_m5');
write_complex_binary(s_m5,'/media/lucasrodes/1DCA91B91A4B64FD/TFG/gnuradio/Flow2/SEND/1 channel/final_rayfad_2');

tic
clear('x_m5','h_m5','s_m5');
x_0 = A(2)*exp(2i*pi*fc/fs*n);
h_0 = h_sigma*(randn(1,S/N)+1i*randn(1,S/N));
h_0 = repelem(h_0,N);
s_0 = h_0.*x_0;
save('Files/Rayleigh Fading/1 Canal/Enviat/signal_0.mat','s_0');
write_complex_binary(s_0,'/media/lucasrodes/1DCA91B91A4B64FD/TFG/gnuradio/Flow2/SEND/1 channel/signal_0');
toc

tic
clear('x_0','h_0','s_0');
x_5 = A(3)*exp(2i*pi*fc/fs*n);
h_5 = h_sigma*(randn(1,S/N)+1i*randn(1,S/N));
h_5 = repelem(h_5,N);
s_5 = h_5.*x_5;
save('Files/Rayleigh Fading/1 Canal/Enviat/signal_5.mat','s_5');
write_complex_binary(s_5,'/media/lucasrodes/1DCA91B91A4B64FD/TFG/gnuradio/Flow2/SEND/1 channel/signal_5');
toc

tic
clear('x_5','h_5','s_5');
x_10 = A(4)*exp(2i*pi*fc/fs*n);
h_10 = h_sigma*(randn(1,S/N)+1i*randn(1,S/N));
h_10 = repelem(h_10,N);
s_10 = h_10.*x_10;
save('Files/Rayleigh Fading/1 Canal/Enviat/signal_10.mat','s_10');
write_complex_binary(s_10,'/media/lucasrodes/1DCA91B91A4B64FD/TFG/gnuradio/Flow2/SEND/1 channel/signal_10');
toc

tic
clear('x_10','h_10','s_10');
x_15 = A(5)*exp(2i*pi*fc/fs*n);
h_15 = h_sigma*(randn(1,S/N)+1i*randn(1,S/N));
h_15 = repelem(h_15,N);
s_15 = h_15.*x_15;
save('Files/Rayleigh Fading/1 Canal/Enviat/signal_15.mat','s_15');
write_complex_binary(s_15,'/media/lucasrodes/1DCA91B91A4B64FD/TFG/gnuradio/Flow2/SEND/1 channel/signal_15');
toc

tic
clear('x_15','h_15','s_15');
x_20 = A(6)*exp(2i*pi*fc/fs*n);
h_20 = h_sigma*(randn(1,S/N)+1i*randn(1,S/N));
h_20 = repelem(h_20,N);
s_20 = h_20.*x_20;
save('Files/Rayleigh Fading/1 Canal/Enviat/signal_20.mat','s_20');
write_complex_binary(s_20,'/media/lucasrodes/1DCA91B91A4B64FD/TFG/gnuradio/Flow2/SEND/1 channel/signal_20');
toc

toc