 clear 
tic
P = 21*10^3; 
Q = 165*10^3; 
a= 3020.5; 
C=8.164; 
w=0.6; 
n=0; 
R0= 10000;
dn = 0.01;
R = gpuArray(-R0:dn:R0);
m = 3.630*10^-4; 
I0=1; 
RO=0; 
L=0;  
I=0;
It = gpuArray(zeros(size(R)));
tn = 0.1;
t0 = 0;
t1 = 3000;
dt = 10;
stepn = 200/dn;
stept = (t1-t0)/dt;
for t = t0:dt:t1 
for n = -100:dn:100    
r = (C.*(-t+a).^w); 
RO = m/(4/3*pi*r^3); 
L = sqrt(complex((4.*(- P.^2.*R.^2 + P.^2.*r.^2 - 2.*P.*Q.*R.*n + 2.*P.*Q.*r.^2 - Q.^2.*n.^2 + Q.^2.*r.^2 + R.^2.*r.^2 - 2.*R.*n.*r.^2 + n.^2.*r.^2))./(P.^2 + 2.*P.*Q + Q.^2 + R.^2 - 2.*R.*n + n.^2))); 

It = It+real(I0*exp(-L.*RO*10)/stepn/stept);
end
toc;
 disp(t);
end 
plot(R,It);

title ('Зависимость суммарной интенсивности I от координаты R') 
xlabel('R, мкм') 
ylabel('Ie, Вт/мкм^2')
toc
