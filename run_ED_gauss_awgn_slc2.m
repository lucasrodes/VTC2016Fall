function [Texp_h0_var,Texp_h0_pdf,Texp_h1_var,Texp_h1_pdf,Tth_var,...
    Tth_h0_pdf,Tth_h1_pdf,Pd_exp,Pfa_exp,Pd_th,Pfa_th] = ...
    run_ED_gauss_awgn_slc2(y1,w1,y2,w2,N)
%RUN_ED_CS_AWGN_SLC2(y1, w1, y2, w2, N) - This function compares theoretical 
%and empirical results of the test statistic probability, i.e. its PDF and 
%the ROC curve, for a Gauss distributed TX signal in an AWGN Channel and 
% using a Square-Law Combiner (SLC) [L=2] at the SU.
%
%  y1, y2: received/sensed signals under H1
%  w1, w2: received/sensed signals under H0
%  N: Number of samples needed to compute an energy value of T(y)

% Default definition of the weighting factors
a1=1;
a2=1;

%% Experimental
% Received power
Prx1 = var(y1);
Prx2 = var(y2);
% Noise standard deviation
sigma_w1 = std(w1); 
sigma_w2 = std(w2); 

% Experimental PDF of the Test Statistic under H0
Y_h0_1 = reshape(w1,N,[]);
T_h0_1 = mean(abs(Y_h0_1).^2,1);
Y_h0_2 = reshape(w2,N,[]);
T_h0_2 = mean(abs(Y_h0_2).^2,1);
T_h0 = a1*T_h0_1+a2*T_h0_2;
[Texp_h0_pdf,Texp_h0_var] = var2pdf(T_h0,200);


% Experimental PDF of the Test Statistic under H1
Y_h1_1 = reshape(y1,N,[]);
T_h1_1 = mean(abs(Y_h1_1).^2,1);
Y_h1_2 = reshape(y2,N,[]);
T_h1_2 = mean(abs(Y_h1_2).^2,1);
T_h1 = a1*T_h1_1+a2*T_h1_2;
[Texp_h1_pdf,Texp_h1_var] = var2pdf(T_h1,200);

disp('Experimental part finished');

%% Theoretical
cont = 0; % A counter variable

% x-axis for the theoretical expression of T(Y) under both hypothesis. Also
% the range of thresholds lambda needed to obtain the ROC curve
Tth_var = min(Texp_h0_var):(max(Texp_h1_var)-min(Texp_h0_var))/1e3:...
    max(Texp_h1_var);

% Theoretical PDF of the Test Statistic under H0
Tth_h0_pdf = pdf('Normal',Tth_var,a1*sigma_w1^2+a2*sigma_w2^2,sqrt(...
    (a1^2*sigma_w1^4+a2^2*sigma_w2^4)/N));

% Theoretical PDF of the Test Statistic under H1
Tth_h1_pdf = pdf('Normal',Tth_var,a1*Prx1+a2*Prx2,sqrt((a1^2*Prx1^2+...
    a2^2*Prx1^2)/N));


% Define the vectors of the Probability of False Alarm (Pfa) and
% Probability of Detection (Pd). Experimental and theoretical.
Pd_exp = zeros(size(Tth_var)); 
Pfa_exp = Pd_exp;
Pd_th = Pd_exp; 
Pfa_th = Pd_exp;

h = waitbar(0,'(SLC AWGN) Please wait...'); % Initiate waitbar
for Tk=Tth_var
    cont = cont + 1; % Update counter
           
    % Pd for threshold Tk (Experimental + Theoretical)
    Pd_exp(cont) = sum(T_h1 > Tk)/length(T_h1);
    Pd_th(cont) = qfunc((Tk-(a1*Prx1+a2*Prx2))./sqrt((a1^2*Prx1^2+...
        a2^2*Prx1^2)/N));
    % Pfa for threshold Tk (Experimental + Theoretical)
    Pfa_exp(cont) = sum(T_h0 > Tk)/length(T_h0);
    Pfa_th(cont) = qfunc((Tk-(a1*sigma_w1^2+a2*sigma_w2^2))/sqrt((a1^2*...
        sigma_w1^4+a2^2*sigma_w2^4)/N));
    
    % Update waitbar
    waitbar(cont/length(Tth_var),h,[num2str(cont/length(Tth_var)*100) ' %']);
end

close(h); % Close waitbar

end

