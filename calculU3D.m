syms x y z t C1 C2 C3 C4 C5 C6 C7 C8 X(x) Y(y) Z(z) W(t) k n m l t coef hknml(t) C bknml(t) b
coefIntegralaSimplificata = zeros(4,4,4,4);
coefFuncW = sym(zeros(1,4));
coefFuncX = sym(zeros(1,4));
coefFuncY = sym(zeros(1,4));
coefFuncZ = sym(zeros(1,4));
bVal = sym(zeros(1,1,1,4));
t0 = 5;
a = 2;
b = 1;
c = 3;
alfa = exp(x) + 2*y + 2*z;
beta = x^2 + 3*y + 2*z;
rezF = beta*t/t0 + alfa*(t0-t)/t0;

X(x) = C1*cos(k*x) + C2*sin(k*x);
Dx = diff(X);
rezX = solve([X(x), X(0)== 1, Dx(0)== 1], [x, C1, C2]);
C1solX = rezX.C1;
C2solX = rezX.C2;
funcX = subs(subs(X(x), C1, C1solX), C2, C2solX);

Y(y) = C3*cos(n*y) + C4*sin(n*y);
Dy = diff(Y);
rezY = solve([Y(y), Y(0)== 1, Dy(0)== 1], [y, C3, C4]);
C3solY = rezY.C3;
C4solY = rezY.C4;
funcY = subs(subs(Y(y), C3, C3solY), C4, C4solY);


Z(z) = C5*cos(m*z) + C6*sin(m*z);
Dz = diff(Z);
rezZ = solve([Z(z), Z(0)== 1, Dz(0)== 1], [z, C5, C6]);
C5solZ = rezZ.C5;
C6solZ = rezZ.C6;
funcZ = subs(subs(Z(z), C5, C5solZ), C6, C6solZ);

W(t) = C7*cos(l*t) + C8*sin(l*t);
Dt = diff(W);
rezW = solve([W(t), W(0)== 1, Dt(0)== 1], [t, C7, C8]);
C7solW = rezW.C7;
C8solW = rezW.C8;
funcW = subs(subs(W(t), C7, C7solW), C8, C8solW);

inmultXYZW = rezF.*funcX.*funcY.*funcZ.*funcW;
calculIntegrala = (int(int(int(int(inmultXYZW,t,0,t0),z,0,c),y,0,b),x,0,a));
integralaSimplificata = simplify(calculIntegrala);

for i = 1:4
    for j = 1:4
        for s = 1:4
            for q = 1:4
                auxIntegralaSimplificata = round(double(simplify(subs(integralaSimplificata,{k,n,m,l}, {i, j, s, q}))),2);
                coefIntegralaSimplificata(i,j,a,b) = auxIntegralaSimplificata(1,1,1,1);
            end
        end
    end
end
display(coefIntegralaSimplificata);

for i = 1:4
     auxFuncX = simplify(subs(funcX,k,i));
     coefFuncX(i) = auxFuncX(i);
end
for i = 1:4
     auxFuncY = simplify(subs(funcY,n,i));
     coefFuncY(i) = auxFuncY(i);
end

for i = 1:4
     auxFuncZ = simplify(subs(funcZ,m,i));
     coefFuncZ(i) = auxFuncZ(i);
end
for i = 1:4
     auxFuncW = simplify(subs(funcW,l,i));
     coefFuncW(i) = auxFuncW(i);
end
for i = 1:4
    for j = 1:4
        for s = 1:4
            for q = 1:4
                hknml(t) = exp(2*int(diff(coefFuncW(1,q))/coefFuncW(1,q)))*(C + int(coefIntegralaSimplificata(1,q)*exp(int(-diff(coefFuncW(1,q))/coefFuncW(1,q)))));
                hknml(t) = subs(hknml(t),C,2);
                display(hknml);
                bknml(t) = int(hknml(t));
                bVal(i,j,s,q) = bknml(t);
            end
        end
    end
end
display(bVal);
