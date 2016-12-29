disp('1.��ʼ���ṹ');
%-----------------------��ʼ���ṹ----------------------------
epsilon0=(1/(36*pi))*1e-9;   %��0=8.854187817*10^-12F/ m  
mu0=4*pi*1e-7;               %��0=4�С�10^-7ţ��/����^2
k0=0;
c=3e+8;
dx=0.05e-6;    %��λm
dy=0.05e-6;
dz=0.05e-6;
dt=0.5*1/(c*sqrt((1/dx)^2+(1/dy)^2+(1/dz)^2));  %��ֵɫɢ��������dt
Npml=12;
nx=190+2*Npml;  %���Ӿ�������x����
ny=190+2*Npml;
nz=21+2*Npml;
nt=2600;   %����ʱ�䲽��
eps_r=ones(nx,ny,nz);
mu_r=ones(nx,ny,nz);
k_r=ones(nx,ny,nz);

%-----------------------���������----------------------------
% LuBiIG����
% epsr_material=4.85;         % �������������ø�������Ľ�糣��
% M_s=1560;                   % ������ı��ʹŻ�ǿ��1560G��1G=10^-3T
% % M_s=1.56;
% gama=1.758e11;              % ����ϵ��rad/(T��s)
% % B_ex=31.7;                  % �ⲿ�ų�T
% B_ex=72;
% wav=150e-6*3.3;
% w=2*pi*2e12;              % �����Ƶ��Ϊ0.5THz
% % w=2*pi*2e12/3.3;
% % w=2*pi*1.3e-6;
% % wav=c*(2*pi)/w;
% % wav=1.3e-6;
% w_m=mu0*gama*M_s*1000;
% % w_m=4*pi*gama*M_s;
% w_ex=gama*B_ex;
% mu_material=(1+w_ex*w_m/(w_ex^2-w^2));
% % mu_material=1;
% k_material=w*w_m/(w_ex^2-w^2);
% % k_material=0;
% mu_all=((w_ex+w_m)^2-w^2)/(w_ex*(w_ex+w_m)-w^2);
% n_eff=sqrt(mu_all*epsr_material);

% w=c/n_eff/wav*2*pi;

% %-----------------------YIG����---------------------------------------
% epsr_material=16;
% M_s=1780;       % 4pMs(G,1 Gs=10^3 A/m)
% dH=0.3;         % ���������߿�Oe��1 Oe=10^3/4pi A/m��
% gama=1.758e11;  % gama=2.8; % ���űȣ�MHz/Gs,�����Ϊ1.758e11��
% % B_ex=31.7;                  % �ⲿ�ų�T
% % B_ex=35;
% % w=2*pi*1e12/1.43;              % �����Ƶ��Ϊ0.5THz
% B_ex=57.414;
% w=2*pi*1.6e12;              % �����Ƶ��Ϊ1.6THz
% % B_ex=30;    % �޸��ⲿ�ų�Ϊ30T�����ڽ�������ʲôͼ��
% % w=2*pi*1.6e12;              % �����Ƶ��Ϊ1.6THz
% % B_ex=53;    % �޸��ⲿ�ų�Ϊ50T�����ڽ�������ʲôͼ��
% % w=2*pi*1.6e12;              % �����Ƶ��Ϊ1.6THz
% % B_ex=36;
% % w=2*pi*1e12;
% % w=2*pi*2e12/3.3;
% % w=2*pi*1.3e-6;
% % wav=c*(2*pi)/w;
% % wav=1.3e-6;
% w_m=mu0*gama*M_s*1000;
% % w_m=4*pi*gama*M_s;
% w_ex=gama*B_ex+1i*gama*mu0*dH*1e3/(4*pi);
% mu_material=(1+w_ex*w_m/(w_ex^2-w^2));
% % mu_material=1;
% k_material=-w*w_m/(w_ex^2-w^2);
% % k_material=0;
% mu_all=((w_ex+w_m)^2-w^2)/(w_ex*(w_ex+w_m)-w^2);
% n_eff=sqrt(mu_all*epsr_material);
% 
% %---------------------------------------------------------------
% N_w=2*pi/(w*dt);
% m_T=50;

% ����Bi:YIG
epsr_material=5.58;
k_material=1.98e-1;
% k_material=0;
mu_material=1;

% wav=3.5e-6;        % m
% wav=1.8e-6;
wav=0.6e-6;

wav_up=3.5e-6;
wav_down=1e-6;

w=2*pi*c/wav;