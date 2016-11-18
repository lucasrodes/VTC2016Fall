S=1e6;

% Define Normal Distributed RV #1
mu_1 = 1;
sigma_1 = 1;
X_1 = sigma_1*(randn(1,S))+mu_1;

% Define Normal Distributed RV #1
mu_2 = 1;
sigma_2 = 1;
X_2 = sigma_2*(randn(1,S))+mu_2;

% Define variable X = max{X_1,X_2} (Empirical and Theoretical)
X = max(X_1,X_2);
[X_pdf, X_var] = var2pdf(X,200); % Obtain empirical CDF (similar to function hist)

X_pdf_th = 1/sqrt(2*pi*sigma_1^2)*exp(-(X_var-mu_1).^2/(2*sigma_1^2)).*...
    (1-qfunc((X_var-mu_2)/sigma_2))+1/sqrt(2*pi*sigma_2^2)...
    *exp(-(X_var-mu_2).^2/(2*sigma_2^2)).*(1-qfunc((X_var-mu_1)/sigma_1));


plot(X_var,X_pdf); hold on;grid on;
plot(X_var,X_pdf_th,'r');  legend('Simulation PDF','Theory PDF'); 
hold off;
