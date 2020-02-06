syms x t C1 C2 C3 C4 X(x) W(t) k n coef hkn(t) C bkn(t)
coefIntegralaSimplificata = zeros(4,4);
coefFuncW = sym(zeros(1,4));
coefFuncX = sym(zeros(1,4));
b = sym(zeros(1,4));
t0 = 5;
a = 2;
alfa = exp(x);
beta = x^2;
rezF = beta*t/t0 + alfa*(t0-t)/t0;

X(x) = C1*cos(k*x) + C2*sin(k*x);
Dx = diff(X);
rezX = solve([X(x), X(0)== 1, Dx(0)== 1], [x, C1, C2]);
C1solX = rezX.C1;
C2solX = rezX.C2;
funcX = subs(subs(X(x), C1, C1solX), C2, C2solX);

W(t) = C3*cos(n*t) + C4*sin(n*t);
Dt = diff(W);
rezW = solve([W(t), W(0)== 1, Dt(0)== 1], [t, C3, C4]);
C3solW = rezW.C3;
C4solW = rezW.C4;
funcW = subs(subs(W(t), C3, C3solW), C4, C4solW);
display(funcW);
inmultXW = rezF.*funcX.*funcW;
calculIntegrala = (int(int(inmultXW,t,0,t0),x,0,a));
integralaSimplificata = simplify(calculIntegrala);

for i = 1:4
    for j = 1:4
        auxIntegralaSimplificata = round(double(simplify(subs(integralaSimplificata,{k,n}, {i, j}))),2);
        coefIntegralaSimplificata(i,j) = auxIntegralaSimplificata(1,1);
    end
end
display(coefIntegralaSimplificata(1,1));
for j = 1:4
     auxFuncW = simplify(subs(funcW,n,j));
     coefFuncW(j) = auxFuncW(j);
end
for i = 1:4
     auxFuncX = simplify(subs(funcX,k,i));
     coefFuncX(i) = auxFuncX(i);
end
sumU2 = 0;
for i = 1:4
    for j = 1:4
        hkn(t) = exp(2*int(diff(coefFuncW(1,j))/coefFuncW(1,j)))*(C + int(coefIntegralaSimplificata(1,j)*exp(int(-diff(coefFuncW(1,j))/coefFuncW(1,j)))));
        hkn(t) = subs(hkn(t),C,2);
        display(hkn);
        bkn(t) = int(hkn(t));
        b(i,j) = bkn(t);
    end
end
display(b);
