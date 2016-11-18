function [pdf,x] = var2pdf(X,n)
%VAR2PDF(X,n) -- Obtains the PDF of a certain vector X

% x: Vector following a certain probability distribution
% pdf: PDF of the vector

    [ht,x] = hist(X,n);
    area = (x(2)-x(1))*sum(ht);
    pdf = ht./area;

end