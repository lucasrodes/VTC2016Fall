function [Texp_h0_var,Texp_h0_pdf,Texp_h1_var,Texp_h1_pdf,Tth_var,...
    Tth_h0_pdf,Tth_h1_pdf,Pd_exp,Pfa_exp,Pd_th,Pfa_th] = ...
    run_ED_gauss_fading_sls2(y1,w1,y2,w2,N)
%RUN_ED_CS_FADING_SLS(y1, w1, y2, w2, N) - This function compares theoretical 
%and empirical results of the test statistic probability, i.e. its PDF and 
%the ROC curve, for a Gauss distributed TX signal in a Rayleigh Fading
%channel and using Square-Law Selector (SLS) at the SU.
%  y1, y2: received/sensed signals under H1
%  w1, w2: received/sensed signals under H0
%  N: Number of samples needed to compute an energy value of T(y)


%% Experimental
sigma_h = 1; % Rayleigh scale parameter
% Received power
Prx1 = var(y1);
Prx2 = var(y2);
% Noise standard deviation
sigma_w1 = std(w1); 
sigma_w2 = std(w2); 
% PU received Power
aP1 = ((Prx1-sigma_w1^2)/(2*sigma_h^2)); % Transmitted power (estimation)
aP2 = ((Prx2-sigma_w2^2)/(2*sigma_h^2)); % Transmitted power (estimation)

% Experimental PDF of the Test Statistic under H0
Y_h0_1 = reshape(w1,N,[]);
T_h0_1 = mean(abs(Y_h0_1).^2,1);
Y_h0_2 = reshape(w2,N,[]);
T_h0_2 = mean(abs(Y_h0_2).^2,1);
T_h0 = max(T_h0_1,T_h0_2);
[Texp_h0_pdf,Texp_h0_var] = var2pdf(T_h0,200);


% Experimental PDF of the Test Statistic under H1
Y_h1_1 = reshape(y1,N,[]);
T_h1_1 = mean(abs(Y_h1_1).^2,1);
Y_h1_2 = reshape(y2,N,[]);
T_h1_2 = mean(abs(Y_h1_2).^2,1);
T_h1 = max(T_h1_1,T_h1_2);
[Texp_h1_pdf,Texp_h1_var] = var2pdf(T_h1,200);

disp('Experimental part finished');

%% Theoretical
Max = 1e2; % Upper bound for integral calculation
cont = 0; % A counter variable

% x-axis for the theoretical expression of T(Y) under both hypothesis. Also
% the range of thresholds lambda needed to obtain the ROC curve
Tth_var = min(Texp_h0_var):(max(Texp_h1_var)-min(Texp_h0_var))/1e3:...
    max(Texp_h1_var);

% Define the vector of the pdf of T(Y)|H1
Tth_h1_pdf = zeros(size(Tth_var));  

% Define the vectors of the Probability of False Alarm (Pfa) and
% Probability of Detection (Pd). Experimental and theoretical.
Pd_exp = zeros(size(Tth_var)); 
Pfa_exp = Pd_exp;
Pd_th = Pd_exp; 
Pfa_th = Pd_exp;

h = waitbar(0,'(SLS Fading) Please wait...'); % Initiate waitbar
for Tk=Tth_var
    cont = cont + 1; % Update counter
    
    % Theoretical PDF of the T(Y)|H1
    f = @(g1,g2) ...
        (1./sqrt(2*pi*1/N*(g1.*aP1+sigma_w1^2).^2)...
        .*exp(-(Tk-(g1.*aP1+sigma_w1^2)).^2./(2*(1/N*(g1.*aP1+...
        sigma_w1^2).^2)))...
        ...
        .*(normcdf((Tk),((g2.*aP2+sigma_w2^2)),sqrt(1/N*(g2.*aP2+...
        sigma_w2^2).^2)))...
        ...
        +1./sqrt(2*pi*1/N*(g2.*aP2+sigma_w2^2).^2)...
        .*exp(-(Tk-(g2.*aP2+sigma_w2^2)).^2./(2*(1/N*(g2.*aP2+...
        sigma_w2^2).^2)))...
        ...
        .*(normcdf((Tk),((g1.*aP1+sigma_w1^2)),sqrt(1/N*(g1.*aP1+...
        sigma_w1^2).^2))))...
        ...
        .*1/(2*sigma_h^2).*exp(-g1./(2*sigma_h^2))...
        .*1/(2*sigma_h^2).*exp(-g2./(2*sigma_h^2));  
    
    Tth_h1_pdf(cont) = integral2(f, 0, Max,0,Max); 
    
    % Pd for threshold Tk (Experimental + Theoretical)
    Pd_exp(cont) = sum(T_h1 > Tk)/length(T_h1);
    ff = @(g1,g2) (1-normcdf(Tk,(g1.*aP1+sigma_w1^2),sqrt(1/N*(g1.*aP1+...
        sigma_w1^2).^2)).*normcdf(Tk,(g2.*aP2+sigma_w2^2),sqrt(1/N*(g2.*aP2+...
        sigma_w2^2).^2))).*1/(2*sigma_h^2).*exp(-g1./(2*sigma_h^2))...
        .*1/(2*sigma_h^2).*exp(-g2./(2*sigma_h^2));
    Pd_th(cont) = integral2(ff, 0, Max,0,Max);
    % Pfa for threshold Tk (Experimental + Theoretical)
    Pfa_exp(cont) = sum(T_h0 > Tk)/length(T_h0);
    Pfa_th(cont) = 1-normcdf(Tk,sigma_w1^2,sigma_w1^2/sqrt(N))*...
        normcdf(Tk,sigma_w2^2,sigma_w2^2/sqrt(N));
    
    % Update waitbar
    waitbar(cont/length(Tth_var),h,[num2str(cont/length(Tth_var)*100) ' %']);
end

% Theoretical PDF of the Test Statistic under H0
mu_h0_1 = sigma_w1^2;
sigma_h0_1 = sigma_w1^2/sqrt(N);
mu_h0_2 = sigma_w2^2;
sigma_h0_2 = sigma_w2^2/sqrt(N);

Tth_h0_pdf = 1/sqrt(2*pi*sigma_h0_1^2)*exp(-(Tth_var-mu_h0_1).^2/...
    (2*sigma_h0_1^2)).*normcdf(Tth_var,mu_h0_2,sigma_h0_2)+1/sqrt(2*...
    pi*sigma_h0_2^2)*exp(-(Tth_var-mu_h0_2).^2/(2*sigma_h0_2^2))...
    .*normcdf(Tth_var,mu_h0_1,sigma_h0_1);

close(h); % Close waitbar

end

