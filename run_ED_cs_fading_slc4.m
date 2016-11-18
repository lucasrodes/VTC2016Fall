function [Texp_h0_var,Texp_h0_pdf,Texp_h1_var,Texp_h1_pdf,Tth_var,...
    Tth_h0_pdf,Tth_h1_pdf,Pd_exp,Pfa_exp,Pd_th,Pfa_th] = ...
    run_ED_cs_fading_slc4(y1,w1,y2,w2,y3,w3,y4,w4,N)
%RUN_ED_CS_FADING_SLC(y1, w1, y2, w2, N) - This function compares theoretical 
%and empirical results of the test statistic probability, i.e. its PDF and 
%the ROC curve, for a Complex Sinusoid TX signal in a Rayleigh Fading
%channel and using Square-Law Combiner (SLC) at the SU.
%  y1, y2: received/sensed signals under H1
%  w1, w2: received/sensed signals under H0
%  N: Number of samples needed to compute an energy value of T(y)

%parpool('local',4);
% Default definition of the weighting factors
a1=1;
a2=1;
a3=1;
a4=0;
%% Experimental
sigma_h = 1; % Rayleigh scale parameter
% Received power
Prx1 = var(y1);
Prx2 = var(y2);
Prx3 = var(y3);
Prx4 = var(y4);% Noise standard deviation
sigma_w1 = std(w1); 
sigma_w2 = std(w2);
sigma_w3 = std(w3); 
sigma_w4 = std(w4); 
% PU reeived Power
aP1 = (Prx1-sigma_w1^2)/(2*sigma_h^2); % Transmitted power (estimation)
aP2 = ((Prx2-sigma_w2^2)/(2*sigma_h^2)); % Transmitted power (estimation)
aP3 = (Prx3-sigma_w3^2)/(2*sigma_h^2); % Transmitted power (estimation)
aP4 = ((Prx4-sigma_w4^2)/(2*sigma_h^2)); % Transmitted power (estimation)

% Experimental PDF of the Test Statistic under H0
Y_h0_1 = reshape(w1,N,[]);
T_h0_1 = mean(abs(Y_h0_1).^2,1);
Y_h0_2 = reshape(w2,N,[]);
T_h0_2 = mean(abs(Y_h0_2).^2,1);
Y_h0_3 = reshape(w3,N,[]);
T_h0_3 = mean(abs(Y_h0_3).^2,1);
Y_h0_4 = reshape(w4,N,[]);
T_h0_4 = mean(abs(Y_h0_4).^2,1);
T_h0 = a1*T_h0_1+a2*T_h0_2+a3*T_h0_3+a4*T_h0_4;
[Texp_h0_pdf,Texp_h0_var] = var2pdf(T_h0,200);


% Experimental PDF of the Test Statistic under H1
Y_h1_1 = reshape(y1,N,[]);
T_h1_1 = mean(abs(Y_h1_1).^2,1);
Y_h1_2 = reshape(y2,N,[]);
T_h1_2 = mean(abs(Y_h1_2).^2,1);
Y_h1_3 = reshape(y3,N,[]);
T_h1_3 = mean(abs(Y_h1_3).^2,1);
Y_h1_4 = reshape(y4,N,[]);
T_h1_4 = mean(abs(Y_h1_4).^2,1);
T_h1 = a1*T_h1_1+a2*T_h1_2+a3*T_h1_3+a4*T_h1_4;
[Texp_h1_pdf,Texp_h1_var] = var2pdf(T_h1,200);

disp('Experimental part finished');

%% Theoretical
Max = 1e2; % Upper bound for integral calculation

% x-axis for the theoretical expression of T(Y) under both hypothesis. Also
% the range of thresholds lambda needed to obtain the ROC curve
delta = (max(Texp_h1_var)-min(Texp_h0_var))/2;
Tth_var = min(Texp_h0_var):delta:max(Texp_h1_var);

% Define the vector of the pdf of T(Y)|H1
Tth_h1_pdf = zeros(size(Tth_var));  

% Define the vectors of the Probability of False Alarm (Pfa) and
% Probability of Detection (Pd). Experimental and theoretical.
Pd_exp = zeros(size(Tth_var)); 
Pfa_exp = Pd_exp;
Pd_th = Pd_exp; 
Pfa_th = Pd_exp;

h = waitbar(0,'(SLC Fading) Please wait...'); % Initiate waitbar
%g3=0; 
%g4=0;
NT = length(Tth_var);
tic
for i = 1:NT
    Tk=min(Texp_h0_var)+i*delta;
    % Theoretical PDF of the T(Y)|H1
    tic
    f = @(g1,g2,g3,g4) 1./(sqrt(2*pi*1/N*(a1^2*(2*g1.*aP1*sigma_w1^2 ...
    +sigma_w1^4)+a2^2*(2*g2.*aP2*sigma_w2^2+sigma_w2^4)...
    +a3^2*(2*g3.*aP3*sigma_w3^2+sigma_w3^4)+a4^2*(2*g4.*aP4*sigma_w4^2+...
    sigma_w4^4))))...
    .*exp(-(Tk-(a1*(g1.*aP1+sigma_w1^2)+a2*(g2.*aP2+sigma_w2^2)...
    +a3*(g3.*aP3+sigma_w3^2)+a4*(g4.*aP4+sigma_w4^2)))...
    .^2./(2*1/N*(a1^2*(2*g1.*aP1*sigma_w1^2+sigma_w1^4)+a2^2*...
    (2*g2.*aP2*sigma_w2^2+sigma_w2^4)+a3^2*(2*g3.*aP3*sigma_w3^2+...
    sigma_w3^4)+a4^2*(2*g4.*aP4*sigma_w4^2+sigma_w4^4)))) ...
    ...
    .*1./(2.*sigma_h^2).*exp(-g1./(2*sigma_h^2))...
    .*1./(2.*sigma_h^2).*exp(-g2./(2*sigma_h^2))...
    .*1./(2.*sigma_h^2).*exp(-g3./(2*sigma_h^2))...
    .*1./(2.*sigma_h^2).*exp(-g4./(2*sigma_h^2)); 
    %Pd_th(i) = integral(@(g4) integral3(@(g1,g2,g3) f, 0, Max, 0, Max, 0, Max), 0, Max,'ArrayValued',true);
    Tth_h1_pdf(i) = integral4(f, 0, Max, 0, Max, 0, Max, 0, Max); 
    %Tth_h1_pdf(i) = integral(@(g4)integral3(@(g1,g2,g3)f(g1,g2,g3,g4),0,Max,0,Max,0,Max),0,Max);
    disp(['T: ' num2str(i)]);
    toc
    % Pd for threshold Tk (Experimental + Theoretical)
    tic
    Pd_exp(i) = sum(T_h1 > Tk)/length(T_h1);
    ff = @(g1,g2,g3,g4) qfunc((Tk-(a1*(g1*aP1+sigma_w1^2)+a2*(g2*aP1+...
    sigma_w2^2)+a3*(g3*aP3+sigma_w3^2)+a4*(g4*aP4+sigma_w4^2)))...
    ./(sqrt((a1^2*(2.*g1*aP1*sigma_w1^2+sigma_w1^4)+...
    a2^2*(2.*g2*aP2*sigma_w2^2+sigma_w2^4)+...
    a3^2*(2.*g3*aP3*sigma_w3^2+sigma_w3^4)+...
    a4^2*(2.*g4*aP4*sigma_w4^2+sigma_w4^4))/N)))...
    .*1./(2.*sigma_h^2).*exp(-g1./(2*sigma_h^2))...
    .*1./(2.*sigma_h^2).*exp(-g2./(2*sigma_h^2))...
    .*1./(2.*sigma_h^2).*exp(-g3./(2*sigma_h^2))...
    .*1./(2.*sigma_h^2).*exp(-g4./(2*sigma_h^2));
    %Pd_th(i) = integral(@(g4) integral3(@(g1,g2,g3) ff, 0, Max, 0, Max, 0, Max), 0, Max,'ArrayValued',true);
    tic
    Pd_th(i) = integral4(ff, 0, Max, 0, Max, 0, Max, 0, Max);
    toc
    tic
    Pd_th(i) = integral(@(g4)integral3(@(g1,g2,g3)ff(g1,g2,g3,g4),0,Max,0,Max,0,Max),0,Max,'ArrayValued',true);
    toc
    
    % Pfa for threshold Tk (Experimental + Theoretical)
    Pfa_exp(i) = sum(T_h0 > Tk)/length(T_h0);
    Pfa_th(i) = qfunc((Tk-(a1*sigma_w1^2+a2*sigma_w2^2+a3*sigma_w3^2 ...
    +a4*sigma_w4^2))/sqrt((a1^2*sigma_w1^4+a2^2*sigma_w2^4 ...
    +a3^2*sigma_w3^4+a4^2*sigma_w4^4)/N));
    disp(['Pd: ' num2str(i)]);
    % Update waitbar
    
    %waitbar(i/length(Tth_var),h,[num2str(i/length(Tth_var)*100) ' %']);
end
toc
disp('fin');
% Theoretical PDF of the Test Statistic under H0
Tth_h0_pdf = pdf('Normal',Tth_var,a1*sigma_w1^2+a2*sigma_w2^2 ...
+a3*sigma_w3^2+a4*sigma_w4^2,sqrt((a1^2*sigma_w1^4+a2^2*sigma_w2^4 ...
+a3^2*sigma_w3^4+a4^2*sigma_w4^4)/N));

close(h); % Close waitbar

end

