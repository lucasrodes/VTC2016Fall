function run_ED(TXsignal, channel, RXtech, N, SNR)
%RUN_ED(TXsignal, channel, RXtech, N, SNR) - 
%  TXsignal: Type of transmitted signal by the PU (Complex sinusoid,
%  Gaussian or PSK)
%  channel: Determines which channel is used (AWGN or Rayleigh Fading)
%  RXtech: 3-sized vector that determines which techniques are used at the
%  SU (Single Antenna, SLC or/and SLS)
%  N: Number of samples needed to compute an energy value of T(y)
%  SNR: SNR under H1 at the SU

% Close windows
try
    close(figure(1));
    close(figure(2));
    close(figure(3));
catch e
end


% Ensure that a RX technique has been selected (Single Antenna, SLC or SLS)
if sum(RXtech)==0
    warndlg('Please select a RX technique');
else
    switch TXsignal
    
    % PU transmitted signal is a complex sinusoid
    case 1
        run_ED_cs(channel, RXtech, N, SNR)
        
    % PU transmitted signal is gaussian distributed
    case 2
        run_ED_gauss(channel, RXtech, N, SNR)
        
    % PU transmitted signal is a psk modulated signal
    case 3
        run_ED_psk(channel, RXtech, N, SNR)
        
    otherwise
        warndl('Unknown error.');
    end
end

