% This script deals as signal generation for Fading channels for a gauss
% distributed PU signal

% Number of Samples
S = 1e7;
% Number of signalsamples with same fading factor
N = 1e3;
% Scale factor of Fading
sigma_h = 1; 
% Transmitted Power (may be changed)
P = 0.05;

% Generate signal
x = sqrt(P)*(randn(1,S)+1i*randn(1,S));
h = sigma_h*(randn(1,S/N)+1i*randn(1,S/N));
h = myrepelem(h,N);
s = h.*x; 

% Save signal 
save('signal_gauss','s');
%write_complex_binary(s,'path/TXsignal_gauss_fading');
