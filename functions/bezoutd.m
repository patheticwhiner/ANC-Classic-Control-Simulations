function [Rp, Sp, nrp, nsp] = bezoutd(A, B, Hs, Hr, P)
%BEZOUTD 求解离散时间Bezout/Diophantine方程
%   通过系数比较法求解Bezout/Diophantine方程 A(z^-1)*Hs(z^-1)*Sp(z^-1) + B(z^-1)*Hr(z^-1)*Rp(z^-1) = P(z^-1)
%   
%   通过求解该方程，可以设计满足极点配置要求的离散时间反馈控制器。
%
% 语法:
%   [Rp, Sp, nrp, nsp] = bezoutd(A, B, Hs, Hr, P)
%
% 输入参数:
%   A   - 系统分母多项式系数向量 [a0 a1 ... aNa] 
%         表示 A(z^-1) = a0 + a1*z^(-1) + ... + aNa*z^(-Na)
%   B   - 系统分子多项式系数向量 [b0 b1 ... bNb]
%         表示 B(z^-1) = b0 + b1*z^(-1) + ... + bNb*z^(-Nb)
%         注意：延迟和离散化延迟需要在B中体现(前导零)
%   Hs  - 控制器分母固定部分系数向量 [hs0 hs1 ... hsNhs]
%         表示 Hs(z^-1) = hs0 + hs1*z^(-1) + ... + hsNhs*z^(-Nhs)
%   Hr  - 控制器分子固定部分系数向量 [hr0 hr1 ... hrNhr]
%         表示 Hr(z^-1) = hr0 + hr1*z^(-1) + ... + hrNhr*z^(-Nhr)
%   P   - 期望闭环特征多项式系数向量 [p0 p1 ... pNp]
%         表示 P(z^-1) = p0 + p1*z^(-1) + ... + pNp*z^(-Np)
%
% 输出参数:
%   Rp  - 求得的控制器分子多项式系数向量 [rp0 rp1 ... rpNrp]
%   Sp  - 求得的控制器分母多项式系数向量 [sp0 sp1 ... spNsp]
%   nrp - Rp多项式的阶数
%   nsp - Sp多项式的阶数
%
% 理论背景:
%   离散时间系统：A(z^-1)y(t) = B(z^-1)u(t)
%   控制律：[Hs(z^-1)*Sp(z^-1)]u(t) = [Hr(z^-1)*Rp(z^-1)][w(t)-y(t)]
%   闭环特征多项式：A(z^-1)*Hs(z^-1)*Sp(z^-1) + B(z^-1)*Hr(z^-1)*Rp(z^-1) = P(z^-1)
%
% 参考文献:
%   I.D. Landau & G. Zito (2006). Digital Control Systems: Design, Identification 
%   and Implementation. Springer.
%
% 作者:
%   J. Langer, I.D. Landau, H. Prochazka
%   7th June 2002
%   修改于 06/08/2007: "convz" 替换为 "conv"

%% 初始化与输入处理
% 计算精度设置
PRECISION = 1e-16;

% 确保所有输入多项式向量为行向量
D = size(A);
if D(1) > 1, A = A'; end

D = size(B);
if D(1) > 1, B = B'; end

D = size(Hs);
if D(1) > 1, Hs = Hs'; end
if D(1) == 0, Hs = 1; end

D = size(Hr);
if D(1) > 1, Hr = Hr'; end
if D(1) == 0, Hr = 1; end

D = size(P);
if D(1) > 1, P = P'; end

%% 计算多项式阶数
na = length(A) - 1;   % A多项式阶数
nb = length(B) - 1;   % B多项式阶数
np = length(P) - 1;   % P多项式阶数
nhs = length(Hs) - 1; % Hs多项式阶数
nhr = length(Hr) - 1; % Hr多项式阶数

%% 计算复合多项式
% 计算A*Hs
if (nhs > 0)
    Ah = conv(A, Hs);
else
    Ah = A * Hs;
end
nah = length(Ah) - 1;

% 计算B*Hr
if (nhr > 0)
    Bh = conv(B, Hr);
else
    Bh = B * Hr;
end
nbh = length(Bh) - 1;

% 检查闭环极点数量
if (np > nah + nbh - 1)
    disp('Bezout错误: 指定的极点过多');
end

%% 处理期望闭环多项式P
% 如果P的阶数不足，补充额外的极点
if (np < nah + nbh - 1)
    % 提取P的现有根
    rootsPdes = roots(P);
    
    % 计算需要添加的额外极点数量
    nextra = nah + nbh - 1 - np;
    
    % 在半径rmin的圆上均匀分布额外极点
    rmin = 1e-16;
    angle = [0:nextra-1]' / nextra * 2 * pi;
    j = sqrt(-1);
    rootsPextra = rmin * exp(j * angle);
    
    % 组合所有根并重构多项式
    P = poly([rootsPdes; rootsPextra]);
    np = nah + nbh - 1;
end
% P, % 显示最终的闭环特征多项式

% 计算控制器多项式阶数
nsp = nbh - 1;
nrp = nah - 1;

%matrix is smaller than vector PD
if (np>nah+nbh-1), 
   disp('模型分母的阶数太低！请将更高阶的多项式添加到 Hs 或 Hr 中。 ');
end;

% ns=nsp+nhs
% nr=nrp+nhr

M=[];
for j=1:nsp+1, 
	V=[];
	if (j>1), V=[V ; zeros(j-1,1)]; end;% zeros in front of Ah
	V=[V ; Ah'];% Ah
	if (j<=nsp), V=[V ; zeros(nsp+1-j,1)]; end;% zeros behind Ah
	if (length(V)~=nah+nbh), disp('bezoutb: error V'); end;
  	M=[M V]; % add one column to M
end;

for j=1:nrp+1, 
	V=[];
	if (j>1), V=[V ; zeros(j-1,1)]; end;
	V=[V ; Bh'];
	if (j<=nrp), V=[V ; zeros(nrp+1-j,1)]; end;
   if (length(V)~=nah+nbh), disp('bezoutb: error V'); end;
  	M=[M V];
end;

D=size(M);
if (D(1)~=nah+nbh), disp('bezoutb: error size M row'); end;
if (D(2)~=nah+nbh), disp('bezoutb: error size M column'); end;

% make P column vector
P=P';


global M1;
M1=M;
      
X= M\P;
% coefficients are real values
X=real(X);

% make X row vector
X=X';
Sp=X(1:nsp+1);
Rp=X(nsp+2:nsp+nrp+2);
