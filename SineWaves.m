clc
close all;
clear all;
figure
%declararea variabilelor

%rezultatul ecuatiei u, care depinde de t
syms u(t)
%variabila utilizata pentru conditiile initiale
syms b 
%k constanta ecuatiei
k = 10;
b = 10;
%ecuatia ce trebuie rezolvata
%d^2u/dt^2 = -k*u, with k = 10
eqn = diff(u,t,2) == -k*u;
%ecuatia diferentiala du/dt
Du = diff(u,t);
%conditiile initiale ale ecuatiei
cond = [u(0)==b, Du(0)==1];
%rezolvarea equatiei 
uSol(t) = dsolve(eqn,cond);
%uSol(t) = C1*cos(t*sqrt(k)+C2*sin(t*sqrt(k)), with
clf;
%reprezentarea grafica a rezultatelor ecuatiei
fplot(uSol);
title('Ecuatia undelor sinusoidale','Color', 'b');
shg;
pause(0.01);
shg;

