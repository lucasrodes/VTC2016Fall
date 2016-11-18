function run_ED_gauss_awgn(RXtech, N, SNR_input)
%RUN_ED_CS_AWGN(RXtech, N, SNR) - This function focuses on a PU gauss 
%distributed signal with an AWGN channel.
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
bb1 = [27 133 165]/255;
bb2 = [133 222 249]/255;
r1 = [130 50 0]/255;
r2 = [230 0 0]/255;
rr1 = [191 8 142]/255;
rr2 = [246 96 206]/255;

% Load directory of Files
files_dir = importdata('Files_Directory.mat');

% Load files
w1 = importdata([files_dir '/Gauss/AWGN/H0_1.mat']);
if SNR_input>=0
    y1 = importdata([files_dir '/Gauss/AWGN/H1_' num2str(SNR_input) ...
        'SNR_1.mat']);
else
    y1 = importdata([files_dir '/Gauss/AWGN/H1_m'...
        num2str(abs(SNR_input)) 'SNR_1.mat']);
end

SNR = 10*log10(var(y1)/var(w1)-1);

if sum(RXtech(2:5))>0 % If any multiple antenna (L=2) technique is selected
    w2 = importdata([files_dir '/Gauss/AWGN/H0_2.mat']);
    if SNR_input>=0
        y2 = importdata([files_dir '/Gauss/AWGN/H1_' num2str(SNR_input)...
            'SNR_2.mat']);
    else
        y2 = importdata([files_dir '/Gauss/AWGN/H1_m' ...
            num2str(abs(SNR_input)) 'SNR_2.mat']);
    end
    
    SNR = 10*log10((var(y1)/var(w1)+var(y2)/var(w2)-2)/2);

    if (RXtech(3)+RXtech(5))>0 % If any multiple antenna (L=4) technique is selected
        w3 = importdata([files_dir '/Gauss/AWGN/H0_3.mat']);
        w4 = importdata([files_dir '/Gauss/AWGN/H0_4.mat']);
        if SNR_input>=0
            y3 = importdata([files_dir '/Gauss/AWGN/H1_' num2str(SNR_input)...
                'SNR_3.mat']);
            y4 = importdata([files_dir '/Gauss/AWGN/H1_' num2str(SNR_input)...
                'SNR_4.mat']);
        else
            y3 = importdata([files_dir '/Gauss/AWGN/H1_m' ...
                num2str(abs(SNR_input)) 'SNR_3.mat']);
            y4 = importdata([files_dir '/Gauss/AWGN/H1_m' ...
                num2str(abs(SNR_input)) 'SNR_4.mat']);
        end
        
        SNR = 10*log10((var(y1)/var(w1)+var(y2)/var(w2)+var(y3)/var(w3)+...
        var(y4)/var(w4)-4)/4);
    
    end
end


%% Three different scenarios: Single Antena, SLC and SLS

% Single Antenna
if RXtech(1) 
    % Run script to obtain all parameters
    [Texp_h0_var,Texp_h0_pdf,Texp_h1_var,Texp_h1_pdf,Tth_var,...
    Tth_h0_pdf,Tth_h1_pdf,Pd_exp,Pfa_exp,Pd_th,Pfa_th] = ...
    run_ED_gauss_awgn_single(y1,w1,N); 

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

% Square-Law Combiner (SLC) [L=2]
if RXtech(2) 
    % Run script to obtain all parameters
    [Texp_h0_var,Texp_h0_pdf,Texp_h1_var,Texp_h1_pdf,Tth_var,...
    Tth_h0_pdf,Tth_h1_pdf,Pd_exp,Pfa_exp,Pd_th,Pfa_th] = ...
    run_ED_gauss_awgn_slc2(y1,w1,y2,w2,N);
    
    % Plot ROC curve
    figure(1);
    plot(Pfa_th,Pd_th,'color',b1,'linewidth',4,'displayname',...
        'Theoretical results (Normal approximation) [AWGN + SLC, L=2]'); 
        hold on;
    plot(Pfa_exp,Pd_exp,'color',b2,'linewidth',2,'displayname',...
        'Experimental results [AWGN + SLC, L=2]'); 
    
    % Plot PDFs of T(Y)
    figure(3)
    plot(Tth_var,Tth_h0_pdf,'color',g1,'linewidth',4); hold on;
    plot(Texp_h0_var,Texp_h0_pdf,'color',g2,'linewidth',2); 
    plot(Tth_var,Tth_h1_pdf,'--','color',g1,'linewidth',4); 
    plot(Texp_h1_var,Texp_h1_pdf,'--','color',g2,'linewidth',2); 
    legend( 'H0 Theoretical (Normal Approximation)','H0 Experimental'...
        ,'H1 Theoretical (Normal Approximation)','H1 Experimental'); 
    title(['AWGN Channel SLC, L=2 (SNR = ' num2str(SNR) 'dB ; N ='...
         num2str(N) ')'],'fontsize',16);
    xlabel('Threshold \lambda','fontsize',16);
    ylabel('PDF of the test statistic T(Y)','fontsize',16);
    grid on; hold off; 

end


% Square-Law Combiner (SLC) [L=4]
if RXtech(3) 
    % Run script to obtain all parameters
    [Texp_h0_var,Texp_h0_pdf,Texp_h1_var,Texp_h1_pdf,Tth_var,...
    Tth_h0_pdf,Tth_h1_pdf,Pd_exp,Pfa_exp,Pd_th,Pfa_th] = ...
    run_ED_gauss_awgn_slc4(y1,w1,y2,w2,y3,w3,y4,w4,N);
    
    % Plot ROC curve
    figure(1);
    plot(Pfa_th,Pd_th,'color',bb1,'linewidth',4,'displayname',...
        'Theoretical results (Normal approximation) [AWGN + SLC, L=4]'); 
        hold on;
    plot(Pfa_exp,Pd_exp,'color',bb2,'linewidth',2,'displayname',...
        'Experimental results [AWGN + SLC, L=4]'); 
    
    % Plot PDFs of T(Y)
    figure(3)
    plot(Tth_var,Tth_h0_pdf,'color',g1,'linewidth',4); hold on;
    plot(Texp_h0_var,Texp_h0_pdf,'color',g2,'linewidth',2); 
    plot(Tth_var,Tth_h1_pdf,'--','color',g1,'linewidth',4); 
    plot(Texp_h1_var,Texp_h1_pdf,'--','color',g2,'linewidth',2); 
    legend( 'H0 Theoretical (Normal Approximation)','H0 Experimental'...
        ,'H1 Theoretical (Normal Approximation)','H1 Experimental'); 
    title(['AWGN Channel SLC, L=4 (SNR = ' num2str(SNR) 'dB ; N ='...
         num2str(N) ')'],'fontsize',16);
    xlabel('Threshold \lambda','fontsize',16);
    ylabel('PDF of the test statistic T(Y)','fontsize',16);
    grid on; hold off; 

end


% Square-Law Selector (SLS) [L=2]
if RXtech(4) 
    % Run script to obtain all parameters
    [Texp_h0_var,Texp_h0_pdf,Texp_h1_var,Texp_h1_pdf,Tth_var,...
    Tth_h0_pdf,Tth_h1_pdf,Pd_exp,Pfa_exp,Pd_th,Pfa_th] = ...
    run_ED_gauss_awgn_sls2(y1,w1,y2,w2,N);
    
    % Plot ROC curve
    figure(1);
    plot(Pfa_th,Pd_th,'color',r1,'linewidth',4,'displayname',...
        'Theoretical results (Normal approximation) [AWGN + SLS, L=2]'); 
        hold on;
    plot(Pfa_exp,Pd_exp,'color',r2,'linewidth',2,'displayname',...
        'Experimental results [AWGN + SLS, L=2]'); 
    
    % Plot PDFs of T(Y)
    figure(3)
    plot(Tth_var,Tth_h0_pdf,'color',g1,'linewidth',4); hold on;
    plot(Texp_h0_var,Texp_h0_pdf,'color',g2,'linewidth',2); 
    plot(Tth_var,Tth_h1_pdf,'--','color',g1,'linewidth',4); 
    plot(Texp_h1_var,Texp_h1_pdf,'--','color',g2,'linewidth',2); 
    legend( 'H0 Theoretical (Normal Approximation)','H0 Experimental'...
        ,'H1 Theoretical (Normal Approximation)','H1 Experimental'); 
    title(['AWGN Channel SLS, L=2 (SNR = ' num2str(SNR) 'dB ; N ='...
         num2str(N) ')'],'fontsize',16);
    xlabel('Threshold \lambda','fontsize',16);
    ylabel('PDF of the test statistic T(Y)','fontsize',16);
    grid on; hold off; 
end


% Square-Law Selector (SLS) [L=2]
if RXtech(5) 
    % Run script to obtain all parameters
    [Texp_h0_var,Texp_h0_pdf,Texp_h1_var,Texp_h1_pdf,Tth_var,...
    Tth_h0_pdf,Tth_h1_pdf,Pd_exp,Pfa_exp,Pd_th,Pfa_th] = ...
    run_ED_gauss_awgn_sls4(y1,w1,y2,w2,y3,w3,y4,w4,N);
    
    % Plot ROC curve
    figure(1);
    plot(Pfa_th,Pd_th,'color',rr1,'linewidth',4,'displayname',...
        'Theoretical results (Normal approximation) [AWGN + SLS, L=4]'); 
        hold on;
    plot(Pfa_exp,Pd_exp,'color',rr2,'linewidth',2,'displayname',...
        'Experimental results [AWGN + SLS, L=4]'); 
    
    % Plot PDFs of T(Y)
    figure(3)
    plot(Tth_var,Tth_h0_pdf,'color',g1,'linewidth',4); hold on;
    plot(Texp_h0_var,Texp_h0_pdf,'color',g2,'linewidth',2); 
    plot(Tth_var,Tth_h1_pdf,'--','color',g1,'linewidth',4); 
    plot(Texp_h1_var,Texp_h1_pdf,'--','color',g2,'linewidth',2); 
    legend( 'H0 Theoretical (Normal Approximation)','H0 Experimental'...
        ,'H1 Theoretical (Normal Approximation)','H1 Experimental'); 
    title(['AWGN Channel SLS, L=4 (SNR = ' num2str(SNR) 'dB ; N ='...
         num2str(N) ')'],'fontsize',16);
    xlabel('Threshold \lambda','fontsize',16);
    ylabel('PDF of the test statistic T(Y)','fontsize',16);
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

