function Dy=GradOp2Dy(f,bc) %Gradient operator with boundary conditions bc
% u(:,:,i_chan)=f 在这里f是矩阵，当i_chan等于1时，就是图像u的第一通道
[r,c]=size(f);
Dy=zeros(r,c,2);
% y1=shift(m,[-1,0],bc)-m; %y(i+1,j)-y(i,j)，这是256*256*2的三阶张量梯度的第一层,是256*256的
% y2=shift(m,[0,-1],bc)-m; %y(i,j+1)-y(i,j)，这是256*256*2的三阶张量梯度的第二层，是256*256的
Gaussian = fspecial('gaussian',[5,5],1.25); %创建二位滤波器(3*3,0.25)
f=conv2(f,Gaussian,'same'); 
y1=shift(f,[-1,0],bc)-f; %y(i+1,j)-y(i,j)，这是256*256*2的三阶张量梯度的第一层,是256*256的
y2=shift(f,[0,-1],bc)-f; %y(i,j+1)-y(i,j)，这是256*256*2的三阶张量梯度的第二层，是256*256的
Tx = 1./ (1 + 13*abs(y1));
Ty = 1./ (1 + 13*abs(y2));
Dy(:,:,1)=Tx;
Dy(:,:,2)=Ty;
end