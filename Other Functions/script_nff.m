y = importdata('h1_f.mat');
w = importdata('H0_1.mat');
N=1e2;

%% Experimental
sigma_h = 1; % Rayleigh fading scale parameter
Prx = var(y); % Received power
sigma_w = std(w); % Noise standard deviation
aP = (Prx-sigma_w^2)/(2*sigma_h^2);

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
fading = 0;

if fading
    Max = 1e2; % Upper bound for integral calculation
    cont = 0; % A counter variable

    % x-axis for the theoretical expression of T(Y) under both hypothesis. Also
    % the range of thresholds lambda needed to obtain the ROC curve
    Tth_var = min(Texp_h0_var):(max(Texp_h1_var)-min(Texp_h0_var))/1e3:...
        max(Texp_h1_var);

    % Define the vector of the pdf
    Tth_h1_pdf = zeros(size(Tth_var));  

    % Define the vectors of the Probability of False Alarm (Pfa) and
    % Probability of Detection (Pd). Experimental and theoretical.
    Pd_exp = zeros(size(Tth_var)); 
    Pfa_exp = Pd_exp;
    Pd_th = Pd_exp; 
    Pfa_th = Pd_exp;

    h = waitbar(0,'(Single Antenna Fading) Please wait...'); % Initiate waitbar
    for Tk=Tth_var
        cont = cont + 1; % Update counter

        % Computing the theoretical expression of T(Y)|H1
        f = @(g) 1./(sqrt(2*pi*1/N*(g.*aP+sigma_w^2).^2))...
            .*exp(-((Tk-(g.*aP+sigma_w^2)).^2./(2*1/N*(g.*aP+sigma_w^2).^2))) ...
         .*1./(2.*sigma_h^2).*exp(-g./(2*sigma_h^2));
        % Analytical pdf for value cont (i.e. Tk)
        Tth_h1_pdf(cont) = integral(f, 0, Max); 

        % Pd for threshold Tk (Experimental + Theoretical)
        Pd_exp(cont) = sum(T_h1 > Tk)/length(T_h1);
        ff = @(g) qfunc((Tk-(g.*aP+sigma_w^2))./(sqrt(1/N*(g.*aP+sigma_w^2).^2))).*1/(2*sigma_h^2).*exp(-g./...
            (2*sigma_h^2));
        Pd_th(cont) = integral(ff, 0, Max);
        % Pfa for threshold Tk (Experimental + Theoretical)
        Pfa_exp(cont) = sum(T_h0 > Tk)/length(T_h0);
        Pfa_th(cont) = qfunc((Tk-sigma_w^2)/(sigma_w^2/sqrt(N)));

        % Update waitbar
        waitbar(cont/length(Tth_var),h,[num2str(cont/length(Tth_var)*100) ' %']);
    end

    % Theoretical PDF of T(Y)|H0
    Tth_h0_pdf = pdf('Normal',Tth_var,sigma_w^2,sqrt(sigma_w^4/N));

    close(h); % Close waitbar
    
else
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
        Pfa_th(cont) = qfunc((Tk-sigma_w^2)./(sigma_w^2/sqrt(N)));

        % Update waitbar
        waitbar(cont/length(Tth_var),h,[num2str(cont/length(Tth_var)*100) ' %']);
    end

    % Theoretical PDF of the Test Statistic under H0
    Tth_h0_pdf = pdf('Normal',Tth_var,sigma_w^2,sigma_w^2/sqrt(N));

    % Theoretical PDF of the Test Statistic under H1
    Tth_h1_pdf = pdf('Normal',Tth_var,Prx,Prx/sqrt(N));

    close(h); % Close waitbar
end

%% Plot
g1 = [51 102 0]/255;
g2 = [0 230 0]/255;
b1 = [0 0 204]/255;
b2 = [102 255 255]/255;
r1 = [130 50 0]/255;
r2 = [230 0 0]/255;

figure(3)
plot(Tth_var,Tth_h0_pdf,'color',g1,'linewidth',4); hold on;
plot(Texp_h0_var,Texp_h0_pdf,'color',g2,'linewidth',2); 
plot(Tth_var,Tth_h1_pdf,'--','color',g1,'linewidth',4); 
plot(Texp_h1_var,Texp_h1_pdf,'--','color',g2,'linewidth',2); 
legend( 'H0 Theoretical (Normal Approximation)','H0 Experimental'...
    ,'H1 Theoretical (Normal Approximation)','H1 Experimental','fontsize',16); 
title(['AWGN Channel (SNR = ' num2str(SNR) 'dB ; N ='...
     num2str(N) ')'],'fontsize',16);
xlabel('Threshold \lambda','fontsize',16);
ylabel('PDF of the test statistic T(Y)','fontsize',16);
grid on; hold off; 