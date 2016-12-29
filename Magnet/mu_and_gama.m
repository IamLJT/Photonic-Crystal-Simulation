% 磁导率张量中mu和k的变化
mu0=4*pi*1e-7;               %μ0=4π×10^-7牛顿/安培^2
epsilon0=(1/(36*pi))*1e-9;   %ε0=8.854187817*10^-12F/ m  
epsilon_r=12.96;
gama=1.758*1e11;            % 单位为rad/(T·s)
Ms=0.156;                   % 单位为T,铁氧体的饱和磁化强度
a=150*1e-6;                 % c/a表示2THz
r=0.35*a;
c=3e8;

w_ex=(100:1000)/100;
w=(1:999)/1000;
w=w';
[W_ex,W]=meshgrid(w_ex,w);
mu=1+W_ex./(W_ex.^2-W.^2);
k=W./(W_ex.^2-W.^2);
Mu=((W_ex+1).^2-W.^2)./(W_ex.*(W_ex+1)-W.^2);
Mu0=(mu.^2-k.^2)./mu;
mesh(W_ex,W,mu);
title('磁导率张量μ');
xlabel('w_e_x(2πc/a)')
ylabel('w(2πc/a)');
zlabel('μ');
figure;
%surf(W_ex,W,mu);
mesh(W_ex,W,k);
title('磁导率张量k');
xlabel('w_e_x(2πc/a)')
ylabel('w(2πc/a)');
zlabel('k');
figure;
%surf(W_ex,W,k);
% figure;
mesh(W_ex,W,Mu);
title('等效磁导率μ_T_E');
xlabel('w_e_x(2πc/a)')
ylabel('w(2πc/a)');
zlabel('μ');
%surf(W_ex,W,Mu);