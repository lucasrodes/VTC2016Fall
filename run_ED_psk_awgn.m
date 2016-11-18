function run_ED_psk_awgn(RXtech, N, SNR)
%RUN_ED_CS_AWGN(RXtech, N, SNR) - This function focuses on a PU PSK 
% modulated signal with an AWGN channel.
%  RXtech: 3-sized vector that determines which techniques are used at the
%  SU (Single Antenna, SLC or/and SLS)
%  N: Number of samples needed to compute an energy value of T(y)
%  SNR: SNR under H1 at the SU

%% Initializations

% Define colors for plot
g1 = [51 102 0]/255;
g2 = [0 230 0]/255;
b1 = [0 0 204]/255;
b2 = [102 255 255]/255;
r1 = [130 50 0]/255;
r2 = [230 0 0]/255;

% Load files
w1 = importdata('Files/PSK/AWGN/H0_1.mat');
if SNR>=0
    y1 = importdata(['Files/PSK/AWGN/H1_' num2str(SNR) ...
        'SNR_1.mat']);
else
    y1 = importdata(['Files/PSK/AWGN/H1_m'...
        num2str(abs(SNR)) 'SNR_1.mat']);
end

if sum(RXtech(2:3))>0 % If any multiple antenna technique is selected
    w2 = importdata('Files/PSK/AWGN/H0_2.mat');
    if SNR>=0
        y2 = importdata(['Files/PSK/AWGN/H1_' num2str(SNR)...
            'SNR_2.mat']);
    else
        y2 = importdata(['Files/PSK/AWGN/H1_m' ...
            num2str(abs(SNR)) 'SNR_2.mat']);
    end
end

%% Three different scenarios: Single Antena, SLC and SLS

% Single Antenna
if RXtech(1) 
    % Run script to obtain all parameters
    [Texp_h0_var,Texp_h0_pdf,Texp_h1_var,Texp_h1_pdf,Tth_var,...
    Tth_h0_pdf,Tth_h1_pdf,Pd_exp,Pfa_exp,Pd_th,Pfa_th] = ...
    run_ED_psk_awgn_single(y1,w1,N); 

    % Plot ROC curve
    figure(1);
    plot(Pfa_th,Pd_th,'color',g1,'linewidth',4,'displayname',...
        'Theoretical results (Normal approximation) [AWGN]');  hold on; 
    plot(Pfa_exp,Pd_exp,'color',g2,'linewidth',2,'displayname',...
        'Experimental results [AWGN]');
    
    % Plot PDFs of T(Y)
    figure(2);
    plot(Tth_var,Tth_h0_pdf,'color',g1,'linewidth',4); hold on;
    plot(Texp_h0_var,Texp_h0_pdf,'color',g2,'linewidth',2); 
    plot(Tth_var,Tth_h1_pdf,'--','color',g1,'linewidth',4); 
    plot(Texp_h1_var,Texp_h1_pdf,'--','color',g2,'linewidth',2);
    legend('H0 Theoretical (Normal Approximation)', ...
        'H0 Experimental','H1 Theoretical (Normal Approximation)',...
        'H1 Experimental'); 
    title(['AWGN Channel (SNR = ' num2str(SNR) 'dB ; N = '...
        num2str(N) ')'],'fontsize',16);
    xlabel('Threshold \lambda','fontsize',16);
    ylabel('PDF of the test statistic T(Y)','fontsize',16);
    grid on; hold off; 

    
end

% Square-Law Combiner (SLC)
if RXtech(2) 
    % Run script to obtain all parameters
    [Texp_h0_var,Texp_h0_pdf,Texp_h1_var,Texp_h1_pdf,Tth_var,...
    Tth_h0_pdf,Tth_h1_pdf,Pd_exp,Pfa_exp,Pd_th,Pfa_th] = ...
    run_ED_psk_awgn_slc(y1,w1,y2,w2,N);
    
    % Plot ROC curve
    figure(1);
    plot(Pfa_th,Pd_th,'color',b1,'linewidth',4,'displayname',...
        'Theoretical results (Normal approximation) [AWGN + SLC]'); 
        hold on;
    plot(Pfa_exp,Pd_exp,'color',b2,'linewidth',2,'displayname',...
        'Experimental results [AWGN + SLC]'); 
    
    % Plot PDFs of T(Y)
    figure(3)
    plot(Tth_var,Tth_h0_pdf,'color',g1,'linewidth',4); hold on;
    plot(Texp_h0_var,Texp_h0_pdf,'color',g2,'linewidth',2); 
    plot(Tth_var,Tth_h1_pdf,'--','color',g1,'linewidth',4); 
    plot(Texp_h1_var,Texp_h1_pdf,'--','color',g2,'linewidth',2); 
    legend( 'H0 Theoretical (Normal Approximation)','H0 Experimental'...
        ,'H1 Theoretical (Normal Approximation)','H1 Experimental'); 
    title(['AWGN Channel SLC (SNR = ' num2str(SNR) 'dB ; N ='...
         num2str(N) ')']);
    xlabel('Threshold \lambda');
    ylabel('PDF of the test statistic T(Y)');
    grid on; hold off; 

     

end

% Square-Law Selector (SLS)
if RXtech(3) 
    % Run script to obtain all parameters
    [Texp_h0_var,Texp_h0_pdf,Texp_h1_var,Texp_h1_pdf,Tth_var,...
    Tth_h0_pdf,Tth_h1_pdf,Pd_exp,Pfa_exp,Pd_th,Pfa_th] = ...
    run_ED_psk_awgn_sls(y1,w1,y2,w2,N);
    
    % Plot ROC curve
    figure(1);
    plot(Pfa_th,Pd_th,'color',r1,'linewidth',4,'displayname',...
        'Theoretical results (Normal approximation) [AWGN + SLS]'); 
        hold on;
    plot(Pfa_exp,Pd_exp,'color',r2,'linewidth',2,'displayname',...
        'Experimental results [AWGN + SLS]'); 
    
    % Plot PDFs of T(Y)
    figure(3)
    plot(Tth_var,Tth_h0_pdf,'color',g1,'linewidth',4); hold on;
    plot(Texp_h0_var,Texp_h0_pdf,'color',g2,'linewidth',2); 
    plot(Tth_var,Tth_h1_pdf,'--','color',g1,'linewidth',4); 
    plot(Texp_h1_var,Texp_h1_pdf,'--','color',g2,'linewidth',2); 
    legend( 'H0 Theoretical (Normal Approximation)','H0 Experimental'...
        ,'H1 Theoretical (Normal Approximation)','H1 Experimental'); 
    title(['AWGN Channel SLS (SNR = ' num2str(SNR) 'dB ; N ='...
         num2str(N) ')']);
    xlabel('Threshold \lambda');
    ylabel('PDF of the test statistic T(Y)');
    grid on; hold off; 
end

% Add title, grid, labels and legend to figure 1 (ROC curve)
figure(1);
xlabel('False Alarm Probability (Pfa)','fontsize',16);
ylabel('Probability of Detection (Pd)','fontsize',16);
title(['ROC Curve (SNR = ' num2str(SNR) 'dB ; N = ' ...
    num2str(N) ')'],'fontsize',16); 
grid on;
hold off;
legend(gca,'show');

end

