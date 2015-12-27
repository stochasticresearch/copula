clear;
clc;

% Test script for estMteDensity
x = normrnd(0,1,1,1000);
[f,xi] = ksdensity(x);
mte_params = estMteDensity(x);

% reconstruct the pdf from the MTE
mte_pdf_estimate = [];
for ii=1:length(mte_params)
    x = mte_params{ii}.xi_subset;
    a = mte_params{ii}.a; b = mte_params{ii}.b;
    c = mte_params{ii}.c; d = mte_params{ii}.d;
    mte_subset_reconstruct = a*exp(b*x)+c*exp(d*x);
    mte_pdf_estimate = [mte_pdf_estimate; mte_subset_reconstruct];
end

plot(xi,f,xi,mte_pdf_estimate);
grid on;
xlabel('x')
legend('KDE', 'MTE')