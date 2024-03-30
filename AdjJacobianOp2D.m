function v = AdjJacobianOp2D(D,w,y,bc)
% Adjoint of the Discrete weighted Jacobian operator(see JacobianOp2D).��ɢ��Ȩ�ſɱ����ӵİ���
% �˴� D Ϊ F���� P �������� F ԭʼ P=zeros([nx ny 2 NGx*NGy*nc])����һ����ά����NG = NGx*NGy;Nc = NcNG/NG;
% D: 4D array with dimensions (Nx,Ny,2,Nc*NG), which for each pixel (i,j)
% contains the generalized 2 x Nc*NG jacobian D(i,j,:,:). The column that
% corresponds to the 2D gradient for the i_sh-th shifting and i_chan
% channel is D(i,j,:, (i_sh-1)*Nc + i_chan)
%
% w: NGx x NGy array with the square root w=sqrt(G) of the convolution
% kernel G (all elements of G must be >=0). NGx , NGy must be odd 
% numbers, since the origin is considered to correspond to the middle
% element of G 
%
% bc: boundary conditions type: 'symmetric' |'circular'|'zero'.
%
% (OUTPUT):
%
% v: Nx x Ny x Nc array (same dims as vectorial image) with the result of
% the adjoint operator on D


if nargin < 3
  bc='symmetric';
end

[Nx,Ny,Nspat_dim,NcNG] = size(D);
[NGx,NGy] = size(w);
NG = NGx*NGy;
Nc = NcNG/NG;%Nc������ͼ���ά��

if Nspat_dim~=2 %ȷ��������ԣ������� A �� B �����ʱ�����Ӧ�����ϵ�Ԫ����Ϊ�߼�ֵ 1 (true)��������Ϊ�߼�ֵ 0
   error('The 3rd dimension of input array D must be 2, since D must have dimensions of the form (Nx,Ny,2,Nc*NG)'); 
end
if floor(Nc) ~= Nc %���������ȡ��
    error('The 4th dimension of input array D must be a multiple of NG, since D must have dimensions of the form (Nx,Ny,2,Nc*NG)');
end
% ��������D�ĵ��ĸ�ά�ȱ�����NG�ı���
Lx = (NGx-1)/2; Ly = (NGy-1)/2;
[shiftsY1,shiftsY2] = ndgrid(-Lx:Lx,-Ly:Ly);

grad_y = zeros(Nx,Ny,2,Nc);
for i_chan=1:Nc
    % gradient of each image channel:
    grad_y(:,:,:,i_chan) = GradOp2Dy(y(:,:,i_chan),bc);
end

F = zeros(Nx,Ny,2,Nc);
for i_sh=1:NG
    
    % adjoint of the weighted stacking of shiftings
    % (single-element indexing of w,shiftsY* corresponds to columnwise
    % scanning of the elements of the 2D arrays, which is what we want)
    % w(i_sh����һ������1��NG������˹������Ԫ�أ�����D�Ĵ�С����F
    F = F + w(i_sh)*shiftAdjST(grad_y.*D(:,:,:,(i_sh-1)*Nc + (1:Nc)), [shiftsY1(i_sh),shiftsY2(i_sh) 0 0],bc);
end
clear D;

v = zeros(Nx,Ny,Nc);
for i_chan=1:Nc    
    v(:,:,i_chan) = AdjGradOp2D(F(:,:,:,i_chan),bc);%F(:,:,:,i_chan)��һ������������������ά��Ϊ2
end



function g=AdjGradOp2D(P,bc) %Adjoint gradient operator (i.e. -div)

P1=P(:,:,1);
P1=shiftAdjST(P1,[-1,0],bc)-P1; % P1(i-1,j)-P1(i,j)
P2=P(:,:,2);
P2=shiftAdjST(P2,[0,-1],bc)-P2; % P2(i,j-1)-P2(i,j)
g=P1+P2;