
%SCRIPT_ED_GAUSS_AWGN_SINGLE(y,w,N) - This script compares theoretical and 
%empirical results of the test statistic probability, i.e. its PDF and the 
%ROC curve, for a Gauss distributed TX signal in an AWGN channel and using 
%a single antenna at the SU.
%  y: received/sensed signal under H1
%  w: received/sensed signal under H0
%  N: Number of samples needed to compute an energy value of T(y)

w = read_complex_binary('path_to_H0_file');
%importdata('Files/Gauss/AWGN/H0_1.mat');
y = read_complex_binary('path_to_H1_file');
%importdata('Files/Gauss/AWGN/H1_0SNR_1.mat');

% Comment if both y and w have already length = 1e7
y(9999873:1e7) = y(end-127:end);
w(9999873:1e7) = w(end-127:end);

N = 100;
SNR = 10*log10(var(y)/var(w)-1);
%% Experimental
Prx = var(y); % Received power
sigma_w = std(w); % Noise standard deviation

% Experimental PDF of the Test Statistic under H0
Y_h0 = reshape(w,N,[]);
T_h0 = mean(abs(Y_h0).^2,1);
[Texp_h0_pdf,Texp_h0_var] = var2pdf(T_h0,200);

% Experimental PDF of the Test Statistic under H1
Y_h1 = reshape(y,N,[]);
T_h1 = mean(abs(Y_h1).^2,1);
[Texp_h1_pdf,Texp_h1_var] = var2pdf(T_h1,200);

disp('Experimental part finished');

%% Theoretical
cont = 0; % A counter variable

% x-axis for the theoretical expression of T(Y) under both hypothesis. Also
% the range of thresholds lambda needed to obtain the ROC curve
Tth_var = min(Texp_h0_var):(max(Texp_h1_var)-min(Texp_h0_var))/1e3:...
    max(Texp_h1_var);

% Define the vectors of the Probability of False Alarm (Pfa) and
% Probability of Detection (Pd). Experimental and theoretical.
Pd_exp = zeros(size(Tth_var)); 
Pfa_exp = Pd_exp;
Pd_th = Pd_exp; 
Pfa_th = Pd_exp;

h = waitbar(0,'(Single Antenna AWGN) Please wait...'); % Initiate waitbar
for Tk=Tth_var
    cont = cont + 1; % Update counter
           
    % Pd for threshold Tk (Experimental + Theoretical)
    Pd_exp(cont) = sum(T_h1 > Tk)/length(T_h1);
    Pd_th(cont) = qfunc((Tk-Prx)./(Prx/sqrt(N)));
    %Pd_th_ex = marcumq(sqrt(2*N*A^2/sigma_w^2),sqrt(2*N*Tk/sigma_w^2),N);
    % Pfa for threshold Tk (Experimental + Theoretical)
    Pfa_exp(cont) = sum(T_h0 > Tk)/length(T_h0);
    Pfa_th(cont) = qfunc((Tk-sigma_w^2)/(sigma_w^2/sqrt(N)));
    %Pfa_th_ex = gammainc(N/sigma_w^2*Tk,N,'upper');
    
    % Update waitbar
    waitbar(cont/length(Tth_var),h,[num2str(cont/length(Tth_var)*100) ' %']);
end

% Theoretical PDF of the Test Statistic under H0
Tth_h0_pdf = pdf('Normal',Tth_var,sigma_w^2,sqrt(sigma_w^4/N));

% Theoretical PDF of the Test Statistic under H1
Tth_h1_pdf = pdf('Normal',Tth_var,Prx,Prx/sqrt(N));

close(h); % Close waitbar

%% PLOTING

% Define colors for plot
g1 = [51 102 0]/255;
g2 = [0 230 0]/255;
b1 = [0 0 204]/255;
b2 = [102 255 255]/255;
r1 = [130 50 0]/255;
r2 = [230 0 0]/255;

% Plot ROC curve
figure(1);
plot(Pfa_th,Pd_th,'color',g1,'linewidth',4,'displayname',...
    'Theoretical results (Normal approximation) [AWGN]');  hold on; 
plot(Pfa_exp,Pd_exp,'color',g2,'linewidth',2,'displayname',...
    'Experimental results [AWGN]');
xlabel('False Alarm Probability (Pfa)','fontsize',16);
ylabel('Probability of Detection (Pd)','fontsize',16);
title(['ROC Curve (SNR = ' num2str(SNR) 'dB ; N = ' ...
    num2str(N) ')'],'fontsize',16); 
grid on;
hold off;

% Plot PDFs of T(Y)
figure(2);
plot(Tth_var,Tth_h0_pdf,'color',g1,'linewidth',4); hold on;
plot(Texp_h0_var,Texp_h0_pdf,'color',g2,'linewidth',2); 
plot(Tth_var,Tth_h1_pdf,'--','color',g1,'linewidth',4); 
plot(Texp_h1_var,Texp_h1_pdf,'--','color',g2,'linewidth',2);
legend('H0 Theoretical (Normal Approximation)', ...
    'H0 Experimental','H1 Theoretical (Normal Approximation)',...
    'H1 Experimental','fontsize',16); 
title(['AWGN Channel (SNR = ' num2str(SNR) 'dB ; N = '...
    num2str(N) ')'],'fontsize',16);
xlabel('Threshold \lambda','fontsize',16);
ylabel('PDF of the test statistic T(Y)','fontsize',16);
grid on; hold off; 

