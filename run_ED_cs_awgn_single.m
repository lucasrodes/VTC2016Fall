function [Texp_h0_var,Texp_h0_pdf,Texp_h1_var,Texp_h1_pdf,Tth_var,...
    Tth_h0_pdf,Tth_h1_pdf,Pd_exp,Pfa_exp,Pd_th,Pfa_th] = ...
    run_ED_cs_awgn_single(y,w,N)
%RUN_ED_CS_AWGN_SINGLE(y,w,N) - This function compares theoretical and 
%empirical results of the test statistic probability, i.e. its PDF and the 
%ROC curve, for a Complex Sinusoid TX signal in an AWGN channel and using a 
%single antenna at the SU.
%  y: received/sensed signal under H1
%  w: received/sensed signal under H0
%  N: Number of samples needed to compute an energy value of T(y)

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
    Pd_th(cont) = qfunc((Tk-Prx)./(sqrt((2*Prx-sigma_w^2)*sigma_w^2/N)));
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
%Tth_h0_pdf_ex = pdf('Gamma',Tth_var,N,(sigma_w^2/N));

% Theoretical PDF of the Test Statistic under H1
Tth_h1_pdf = pdf('Normal',Tth_var,Prx,sqrt((2*Prx*sigma_w^2-sigma_w^4)/N));
%Tth_h1_pdf_ex = 2*N/sigma_w^2*pdf('Noncentral Chi-square',2*N/sigma_w^2*Tth_var,2*N,2*N*A^2/sigma_w^2);

close(h); % Close waitbar

end



