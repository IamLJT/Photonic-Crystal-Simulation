function e=hole_matrix(epsilon,i,j,k,L,epsilon0,nxp_s,nxp_e,nyp_s,nyp_e)
%	epsilon   ��糣������
%   i         Բ����x����
%   j         Բ����y����
%   k         Բ����z����
%   L         Բ�İ뾶
%   epsilon0  ������糣��
for a=(i-L):(i+L)   %   ʵ�ʰ뾶ΪL+1��������
    for b=(j-L):(j+L)
        if (a-i)^2+(b-j)^2<L^2+0.5&&(a>=nxp_s&&a<=nxp_e&&b>=nyp_s&&b<=nyp_e)
            epsilon(a,b,k)=epsilon0;
        end
    end
end
e=epsilon;