function D = JacobianOp2D(u,w,y,bc)
% Operator of the discrete weighted Jacobian (and related functionals) 
%
% u: Nx x Ny x Nc array of a vector-valued image with Nc-channels, defined
% on a Nx x Ny pixel grid.256*256*3
% 此处 u=project(K)而 K=y-lambda*AdjJacobianOp2D(F,w,bc)
% w: NGx x NGy array with the square root w=sqrt(G) of the convolution
% kernel G (all elements of G must be >=0). NGx , NGy must be odd 
% numbers, since the origin is considered to correspond to the middle
% element of G. 
%
% bc: boundary condition type: 'symmetric' |'circular'|'zero'.
%
% (OUTPUT):
% D: 4D array with dimensions (Nx,Ny,2,Nc*NG)即256*256*2*?, which for each pixel (i,j)
% contains the generalized 2 x Nc*NG Jacobian D(i,j,:,:). The column that
% corresponds to the 2D gradient for the i_sh-th shifting and i_chan
% channel is D(i,j,:, (i_sh-1)*Nc + i_chan)
%一个张量求完雅可比矩阵

if nargin < 3
  bc='symmetric';
end

[Nx,Ny,Nc] = size(u);
% P = Nx*Ny;

[NGx,NGy] = size(w);
NG = NGx*NGy;%NG等于9

if ~all(mod([NGx,NGy], 2)) % if not all [NGx,NGy] are odd numbers
    error('The dimensions of the kernel G must both be odd numbers');
end

Lx = (NGx-1)/2; Ly = (NGy-1)/2;
[shiftsY1,shiftsY2] = ndgrid(-Lx:Lx,-Ly:Ly);


grad_c = zeros(Nx,Ny,2,Nc);
for i_chan=1:Nc
    % gradient of each image channel:
    grad_c(:,:,:,i_chan) = GradOp2D(u(:,:,i_chan),bc);
end
% grad_c的每条通道都是256*256*2的三阶张量，GradOp2D(u(:,:,i_chan),bc)是一个256*256*2的张量
% 当i_chan等于1时有(Nx,Ny,1,1),(Nx,Ny,2,1)这两个数组
% clear u;

grad_y = zeros(Nx,Ny,2,Nc);
for i_chan=1:Nc
    % gradient of each image channel:
    grad_y(:,:,:,i_chan) = GradOp2Dy(y(:,:,i_chan),bc);
end

D = zeros(Nx,Ny,2,Nc*NG);
for i_sh=1:NG
    
    % shift wrt spatial dimensions X and Y:
    % (single-element indexing of w,shiftsY* corresponds to columnwise
    % scanning of the elements of the 2D arrays, which is what we want)
    D(:,:,:, (i_sh-1)*Nc + (1:Nc)) = w(i_sh)*shift(grad_y.*grad_c, [shiftsY1(i_sh),shiftsY2(i_sh) 0 0],bc);
end


