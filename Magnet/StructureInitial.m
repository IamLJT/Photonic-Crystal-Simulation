disp('1.初始化结构');
%-----------------------初始化结构----------------------------
epsilon0=(1/(36*pi))*1e-9;   %ε0=8.854187817*10^-12F/ m  
mu0=4*pi*1e-7;               %μ0=4π×10^-7牛顿/安培^2
k0=0;
c=3e+8;
dx=0.05e-6;    %单位m
dy=0.05e-6;
dz=0.05e-6;
dt=0.5*1/(c*sqrt((1/dx)^2+(1/dy)^2+(1/dz)^2));  %数值色散条件计算dt
Npml=12;
nx=190+2*Npml;  %光子晶体区域x方向
ny=190+2*Npml;
nz=21+2*Npml;
nt=2600;   %运行时间步长
eps_r=ones(nx,ny,nz);
mu_r=ones(nx,ny,nz);
k_r=ones(nx,ny,nz);

%-----------------------铁氧体参数----------------------------
% LuBiIG材料
% epsr_material=4.85;         % 待查阅资料设置该铁氧体的介电常数
% M_s=1560;                   % 铁氧体的饱和磁化强度1560G，1G=10^-3T
% % M_s=1.56;
% gama=1.758e11;              % 旋磁系数rad/(T·s)
% % B_ex=31.7;                  % 外部磁场T
% B_ex=72;
% wav=150e-6*3.3;
% w=2*pi*2e12;              % 入射光频率为0.5THz
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

% %-----------------------YIG材料---------------------------------------
% epsr_material=16;
% M_s=1780;       % 4pMs(G,1 Gs=10^3 A/m)
% dH=0.3;         % 铁磁吸收线宽（Oe，1 Oe=10^3/4pi A/m）
% gama=1.758e11;  % gama=2.8; % 旋磁比（MHz/Gs,换算后即为1.758e11）
% % B_ex=31.7;                  % 外部磁场T
% % B_ex=35;
% % w=2*pi*1e12/1.43;              % 入射光频率为0.5THz
% B_ex=57.414;
% w=2*pi*1.6e12;              % 入射光频率为1.6THz
% % B_ex=30;    % 修改外部磁场为30T，不在禁带内是什么图像？
% % w=2*pi*1.6e12;              % 入射光频率为1.6THz
% % B_ex=53;    % 修改外部磁场为50T，不在禁带内是什么图像？
% % w=2*pi*1.6e12;              % 入射光频率为1.6THz
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

% 材料Bi:YIG
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