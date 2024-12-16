w_std = sqrt(10);
v_std = sqrt(1);
x=randn(1);
y=[];
N_steps=15;
alp=1;
bet=5;

for k=1:N_steps
    w=randn(1)*w_std;
    x(k+1)=0.5*x(k)+bet*x(k)/(1+x(k)^2)+8*cos(1.2*k)+w;
    v=randn(1)*v_std;
    y(k)=alp*x(k)+x(k)^2/20+v;
end

x=x(1:end-1);