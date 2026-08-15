function UXY=testU(M1,N1,v)
N0=512;
       xw=gauss(2*N0+1,r_jacobi(2*N0+1,0,0)); 
       yw=gauss(2*N0+1, r_jacobi(2*N0+1,0,0));
        xi = xw(:,1)'; wr = xw(:,2); 
        yi=yw(:,1)';wt=yw(:,2);
        xa = 0; xb = 1; h = xb - xa; 
        x = 1/2*(xb-xa)*xi + 1/2*(xb+xa);   
        dxdxi = h/2; dxidx = 2/h; 
        ya=-pi/4;yb=3*pi/2;h2=yb-ya;
        y=1/2*(yb-ya)*yi+1/2*(yb+ya);
        dydyi=h2/2;dyidy=2/h2;
        % ==================== 2. 几何定义 ====================
        % 边界形状函数 R(θ)
        mtheta = (cos(y).^4+sin(y).^4).^(1./4);
        mtheta1 = -(2.^(1/2).*sin(4.*y))./(2.*(cos(4.*y) + 3).^(3/4));           % 一阶导数   
        mtheta2 = -(2.^(1/2).*(12.*cos(4.*y) + cos(4.*y).^2 + 3))./(2.*(cos(4.*y) + 3).^(7/4));          % 二阶导数
        p1= 1- mtheta2./mtheta + 2.*mtheta1.^2./mtheta.^2 ;
        p2= 1+ mtheta1.^2./mtheta.^2 ;
        p3= 2.*mtheta1./mtheta;
        p4=1;
        p5= mtheta.^2;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 真解信息
fxy=@(x,y)1;
[xn,ym]=meshgrid(x,y);
FXY=fxy(xn,ym); 
      

% ==================== 3. 基函数构造 ====================
        % Legendre多项式计算
        Le = zeros(N1+1, 2*N0+1); Le1 = Le; Le2 = Le;
        Le(1,:) = 1; Le(2,:) = xi; Le1(2,:) = 1;
         for n = 1:N1
            % Legendre多项式递推
            Le(n+2,:) = ((2*n+1)*xi.*Le(n+1,:) - n*Le(n,:))/(n+1); 
            % 一阶导数
            Le1(n+2,:) = (2*n+1)/(n+1)*Le(n+1,:) + ((2*n+1)*xi.*Le1(n+1,:) - n*Le1(n,:))/(n+1);
            % 二阶导数  
            Le2(n+2,:) = 2*(2*n+1)/(n+1)*Le1(n+1,:) + ((2*n+1)*xi.*Le2(n+1,:) - n*Le2(n,:))/(n+1);
         end      
        phi_ko =  Le(1:N1-1,:) - Le(3:N1+1,:);
        phi_ko1 = Le1(1:N1-1,:) - Le1(3:N1+1,:);
        phi_ko2 = Le2(1:N1-1,:) - Le2(3:N1+1,:);
phi_k=zeros(N1-3,2*N0+1);phi_k1=phi_k;phi_k2=phi_k;
       for k = 1:N1-3
                phi_k(k,:) = phi_ko(k+2,:) -  (2*(k+1)+3)/(2*(k-1)+3) *phi_ko(k,:);
                phi_k1(k,:) =phi_ko1(k+2,:) - (2*(k+1)+3)/(2*(k-1)+3) *phi_ko1(k,:);
                phi_k2(k,:) =phi_ko2(k+2,:) - (2*(k+1)+3)/(2*(k-1)+3) *phi_ko2(k,:);
        end
le = zeros(M1+1, 2*N0+1); 
le1 = le; 
le2 = le;
le(1,:) = 1; 
le(2,:) = yi; 
le1(2,:) = 1;

for n = 1:M1
    le(n+2,:) = ((2*n+1)*yi.*le(n+1,:) - n*le(n,:))/(n+1); 
    le1(n+2,:) = (2*n+1)/(n+1)*le(n+1,:) + ((2*n+1)*yi.*le1(n+1,:) - n*le1(n,:))/(n+1);
    le2(n+2,:) = 2*(2*n+1)/(n+1)*le1(n+1,:) + ((2*n+1)*yi.*le2(n+1,:) - n*le2(n,:))/(n+1);
end
       eta_ko = le(1:M1-1,:) - le(3:M1+1,:);
       eta_ko1 =le1(1:M1-1,:) - le1(3:M1+1,:);
       eta_ko2 =le2(1:M1-1,:) - le2(3:M1+1,:);
eta_k=zeros(M1-3,2*N0+1);eta_k1=eta_k;eta_k2=eta_k;
        for k = 1:M1-3
                eta_k(k,:) = eta_ko(k+2,:) -   (2*(k+1)+3)/(2*(k-1)+3)*eta_ko(k,:);
                eta_k1(k,:) = eta_ko1(k+2,:) - (2*(k+1)+3)/(2*(k-1)+3) *eta_ko1(k,:);
                eta_k2(k,:) = eta_ko2(k+2,:) - (2*(k+1)+3)/(2*(k-1)+3) *eta_ko2(k,:);
        end
        % ==================== 4. 径向矩阵组装 ====================
        % 初始化径向矩阵
        A1 = zeros(N1-3, N1-3); A2 = A1; A3 = A1; A4 = A1; A5 = A1;
        A6 = A1; A7 = A1; A8 = A1; A9 = A1; A10=A1;A11=A1;A12=A1;A13=A1;A14=A1;A15=A1;A16=A1;
        
        for  k= 1:N1-3
            for j = 1:N1- 3
                % 各项积分计算
                A1(k,j) = ((xi+1)./2).^(-1) .* phi_k1(j,:) .* phi_k1(k,:) * wr *dxidx^2* dxdxi;
                A2(k,j) =  phi_k1(j,:) .* phi_k2(k,:) * wr *dxidx^3*dxdxi;
                A3(k,j) = ((xi+1)/2).^(-1) .* phi_k1(j,:) .* phi_k1(k,:) * wr*dxidx^2 * dxdxi;
                A4(k,j) = ((xi+1)./2).^(-2) .* phi_k1(j,:) .* phi_k(k,:) * wr *dxidx* dxdxi;
                A5(k,j) =   phi_k2(j,:) .* phi_k1(k,:) * wr *dxidx^3*  dxdxi;
                A6(k,j) =  ((xi+1)./2) .* phi_k2(j,:) .* phi_k2(k,:) * wr *dxidx^4*  dxdxi;
                A7(k,j) =  phi_k2(j,:) .* phi_k1(k,:) * wr *dxidx^3*  dxdxi;
                A8(k,j) =  ((xi+1)./2).^(-1) .* phi_k2(j,:) .* phi_k(k,:) * wr *dxidx^2*  dxdxi;
                A9(k,j) =  ((xi+1)./2).^(-1) .* phi_k1(j,:) .* phi_k1(k,:) * wr *dxidx^2*  dxdxi;
                A10(k,j)=  phi_k1(j,:) .* phi_k2(k,:) * wr *dxidx^3*  dxdxi;
                A11(k,j)=  ((xi+1)./2).^(-1) .* phi_k1(j,:) .* phi_k1(k,:) * wr *dxidx^2*  dxdxi;
                A12(k,j)=  ((xi+1)./2).^(-2) .* phi_k1(j,:) .* phi_k(k,:) * wr *dxidx*  dxdxi;
                A13(k,j)= ((xi+1)./2).^(-2) .* phi_k(j,:) .* phi_k1(k,:) * wr *dxidx*  dxdxi;
                A14(k,j)= ((xi+1)./2).^(-1) .* phi_k(j,:) .* phi_k2(k,:) * wr *dxidx^2*  dxdxi;
                A15(k,j)= ((xi+1)./2).^(-2) .* phi_k(j,:) .* phi_k1(k,:) * wr *dxidx*  dxdxi;
                A16(k,j)= ((xi+1)./2).^(-3) .* phi_k(j,:) .* phi_k(k,:) * wr *  dxdxi;
            end
        end
        % ==================== 5. 角向矩阵组装 ====================
        % 初始化角向矩阵
D1=zeros(M1-3,M1-3);D2=D1;D3=D1;D4=D1;D5=D1;D6=D1;D7=D1;D8=D1;D9=D1;D10=D1;D11=D1;D12=D1;D13=D1;D14=D1;D15=D1;D16=D1;
for m = 1:M1-3
    for n = 1:M1-3
        % 第一项
                D1(m,n) =   (p1.^2)./p5.*eta_k(n,:) .* eta_k(m,:) * wt * dydyi;
                D2(m,n) =   p1.*p2./p5.*eta_k(n,:) .* eta_k(m,:) * wt *dydyi;
                D3(m,n) =   p1.*p3./p5.*eta_k(n,:) .* eta_k1(m,:) * wt*dyidy * dydyi;
                D4(m,n) =   p1.*p4./p5.*eta_k(n,:) .* eta_k2(m,:) * wt*dyidy^2 * dydyi;
                D5(m,n) =   p2.*p1./p5.*eta_k(n,:) .* eta_k(m,:) * wt *  dydyi;
                D6(m,n) =   p2.*p2./p5.*eta_k(n,:) .* eta_k(m,:) * wt *  dydyi;
                D7(m,n) =   p2.*p3./p5.*eta_k(n,:) .* eta_k1(m,:) * wt *dyidy*  dydyi;
                D8(m,n) =   p2.*p4./p5.*eta_k(n,:) .* eta_k2(m,:) * wt *dyidy^2*  dydyi;
                D9(m,n) =   p3.*p1./p5.*eta_k1(n,:) .* eta_k(m,:) * wt *dyidy*  dydyi;
                D10(m,n)=   p3.*p2./p5.*eta_k1(n,:) .* eta_k(m,:) * wt *dyidy*  dydyi;
                D11(m,n)=   (p3.^2)./p5.*eta_k1(n,:) .* eta_k1(m,:) * wt *dyidy^2*  dydyi;
                D12(m,n)=   p3.*p4./p5.*eta_k1(n,:) .* eta_k2(m,:) * wt *dyidy^3*  dydyi;
                D13(m,n)=   p4.*p1./p5.*eta_k2(n,:) .* eta_k(m,:) * wt *dyidy^2*  dydyi;
                D14(m,n)=   p4.*p2./p5.*eta_k2(n,:) .* eta_k(m,:) * wt *dyidy^2*  dydyi;
                D15(m,n)=   p4.*p3./p5.*eta_k2(n,:) .* eta_k1(m,:) * wt *dyidy^3*  dydyi;
                D16(m,n)=   (p4^2)./p5.*eta_k2(n,:) .* eta_k2(m,:) * wt *dyidy^4*  dydyi;
    end
end

      % 投影到基函数上
     F1 = zeros(M1-3, N1-3); 
        for m = 1:M1-3
            for n = 1:N1-3
                 F1(m,n)=(eta_k(m,:).*wt')*dydyi*FXY*((phi_k(n,:))'.*wr)*dxdxi;
            end
        end 
        
        % ==================== 7. 全局矩阵组装 ====================
        % 使用Kronecker积组装全局矩阵
% 组装全局矩阵
W1 =v.*(kron(A1,D1)+kron(A2,D2)-kron(A3,D3)+kron(A4,D4)+kron(A5,D5)+kron(A6,D6)-kron(A7,D7)+kron(A8,D8)-kron(A9,D9)-kron(A10,D10)+kron(A11,D11)-kron(A12,D12)+kron(A13,D13)+kron(A14,D14)-kron(A15,D15)+kron(A16,D16));
% 右端向量
       F = reshape(F1, [], 1);
  % 原来的
XXX = W1 \ F;

        % ==================== 9. 重构解和后处理 ===================
        XXX1=reshape(XXX,M1-3,N1-3);
        UXY=eta_k'*XXX1*phi_k;