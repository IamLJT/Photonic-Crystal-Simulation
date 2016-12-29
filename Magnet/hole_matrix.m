function e=hole_matrix(epsilon,i,j,k,L,epsilon0,nxp_s,nxp_e,nyp_s,nyp_e)
%	epsilon   介电常数矩阵
%   i         圆中心x坐标
%   j         圆中心y坐标
%   k         圆中心z坐标
%   L         圆的半径
%   epsilon0  空气介电常数
for a=(i-L):(i+L)   %   实际半径为L+1个像素了
    for b=(j-L):(j+L)
        if (a-i)^2+(b-j)^2<L^2+0.5&&(a>=nxp_s&&a<=nxp_e&&b>=nyp_s&&b<=nyp_e)
            epsilon(a,b,k)=epsilon0;
        end
    end
end
e=epsilon;