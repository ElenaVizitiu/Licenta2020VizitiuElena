% lungimea corzii
L = 1;  
% viteza undelor
c = 1;  
% frecventa undei de intrare
f = 5; 
pi = 3.14;
% timpul total de simulare
time = 4; 
% pas în pozitie
dx = 0.005; 
% pasul de timp
dt = 0.005;  
% numarul total de coloane
nx = ceil(L/dx)+1;  
% numarul total de linii
nt = ceil(time/dt)+1;  
% initializarea matricii u cu zerouri
u = zeros(nt,nx);   
r = (c*dt/dx)^2;

% implementarea conditiilor la limita pentru x = 1
for i = 2:nt 
  u(i,1) = 3*sin(2*pi*f*(i-1)*dt)*exp((i-1)*dt);
end 
 
% implementarea ec. 12
for j = 2:nx-1 
   %solutie pentru prima etapa
  u(2,j) = u(1,j)+1/2*r*(u(1,j+1) - 2*u(1,j)+ u(1,j-1)); 
end 

% popularea solutiilor matricii(implementarea ec. 9)
for i = 3:nt-1 
    for j= 2:nx-1 
     u(i,j) = 2*u(i-1,j)-u(i-2,j)+r*((u(i-1,j+1) - 2*u(i-1,j)+ u(i-1,j-1)));  
    end
end

% reprezentarea grafica
[nt,nx] = size(u);
figure
hold on;
for i = 3 : nt
    plot(u(i,:));
    title('Ecuatia de unda unidimensionala','Color', 'b');
    axis([0 nx -5 5]);
    pause(0.0100);
    hold off;
end
