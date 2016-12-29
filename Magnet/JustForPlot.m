c=3e8;
mu0=4*pi*1e-7;               %��0=4�С�10^-7ţ��/����^2
B_ex=linspace(55,57,100);
M_s=1560;                   % ������ı��ʹŻ�ǿ��1560G��1G=10^-3T
gama=1.758e11;              % ����ϵ��rad/(T��s)
w=2*pi*1.6e12;
% w=2*pi*5e14;
% w=2*pi*2e12/5.7143;   
wav=c*(2*pi)/w;
% wav=1.3e-6;
w_m=mu0*gama*M_s*1000;
w_ex=gama*B_ex+1j*gama*300*1e-7/2;
mu_material=(1+w_ex*w_m./(w_ex.^2-w^2));
% mu_material=1;
k_material=w*w_m./(w_ex.^2-w^2);
% k_material=0;
mu_all=((w_ex+w_m).^2-w^2)./(w_ex.*(w_ex+w_m)-w^2);
mu_all_1=(mu_material.^2-k_material.^2)./mu_material;
figure(1);
plot(B_ex,real(mu_material),'r');
hold on;
plot(B_ex,real(k_material),'k');
hold on;
plot(B_ex,real(mu_all),'-.');
hold on;
plot(B_ex,imag(mu_all),'--');
% figure(2);
% plot(B_ex,real(mu_all_1));
% hold on;
% plot(B_ex,imag(mu_all_1),'--');