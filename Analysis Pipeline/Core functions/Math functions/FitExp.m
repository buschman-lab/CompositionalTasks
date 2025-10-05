function [ft,gof]=FitExp(x,y,StartPoint) % fits an exponential function
%ft,gof] = fit(x,y,['exp' num2str(ExpNTerms)] );


%ft = fittype( 'a+b*exp(c*(x-x0))', 'independent', {'x'}, 'dependent', 'y','coefficients',{'a','b','c','x0'});
%ft = fittype( 'ExpFunc(x,a,b,c,x0)', 'independent', {'x'}, 'dependent', 'y','coefficients',{'a','b','c','x0'});
g = fittype('a-b*exp(-c*x)');

if nargin<3
    [ft,gof] = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
else
    [ft,gof] = fit( x, y, g ,'StartPoint', StartPoint );
end