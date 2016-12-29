% �ŵ���������mu��k�ı仯
mu0=4*pi*1e-7;               %��0=4�С�10^-7ţ��/����^2
epsilon0=(1/(36*pi))*1e-9;   %��0=8.854187817*10^-12F/ m  
epsilon_r=12.96;
gama=1.758*1e11;            % ��λΪrad/(T��s)
Ms=0.156;                   % ��λΪT,������ı��ʹŻ�ǿ��
a=150*1e-6;                 % c/a��ʾ2THz
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
title('�ŵ���������');
xlabel('w_e_x(2��c/a)')
ylabel('w(2��c/a)');
zlabel('��');
figure;
%surf(W_ex,W,mu);
mesh(W_ex,W,k);
title('�ŵ�������k');
xlabel('w_e_x(2��c/a)')
ylabel('w(2��c/a)');
zlabel('k');
figure;
%surf(W_ex,W,k);
% figure;
mesh(W_ex,W,Mu);
title('��Ч�ŵ��ʦ�_T_E');
xlabel('w_e_x(2��c/a)')
ylabel('w(2��c/a)');
zlabel('��');
%surf(W_ex,W,Mu);