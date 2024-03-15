function [InterSect,Line1,Line2] = FindLineInterSection(X,F1,F2)
%FINDLINEINTERSECTION finds intersection between two curve. 
ShowPlot=0;
% X is x axis 
% Y1 values of curve 1 and Y2 values of curve 2
P_B=find(diff(F1>=F2))+1;
P_S=find(diff(F1<F2));

X1=X(P_B);X2=X(P_S);
[Line1,m1,a1]=FitLine(X1,F1(P_B),X2,F1(P_S));
[Line2,m2,a2]=FitLine(X1,F2(P_B),X2,F2(P_S));
% 
% IntFunc=@(x) Line1(x)-Line2(x);
% InterSect=fzero(IntFunc,X1);
InterSect=(a2-a1)/(m1-m2);
if ShowPlot
    figure
    hold on
    plot(X,F1,X,F2)
    plot(log(InterSect),Line1(InterSect),'r*')
   % plot(InterSect,Line2(InterSect),'g*')
     set(gca, 'xscale', 'log', 'ytick', [0 0.25 0.5 0.75 1])

end

end
function [L,m,a]=FitLine(X1,Y1,X2,Y2)
% fits a line to two points Y=mx+a
m=(Y2-Y1)/(X2-X1);
a=Y1-m*X1;

L=@(x) m*x+a;


end
