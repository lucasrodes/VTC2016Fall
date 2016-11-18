function run_ED_psk(channel, RXtech, N, SNR)
%RUN_ED_CS(channel, RXtech, N, SNR) - this function focuses on a PU PSK 
% modulated signal.
%  channel: Determines which channel is used (AWGN or Rayleigh Fading)
%  RXtech: 3-sized vector that determines which techniques are used at the
%  SU (Single Antenna, SLC or/and SLS)
%  N: Number of samples needed to compute an energy value of T(y)
%  SNR: SNR under H1 at the SU

switch channel
    
    % AWGN Channel
    case 1
        run_ED_psk_awgn(RXtech,N,SNR);
    
    % Rayleigh Fading Channel
    case 2
        run_ED_psk_fading(RXtech,N,SNR);
        
    otherwise
        warndlg('Unknown error');
end

end

