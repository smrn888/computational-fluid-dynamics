clc
clear all
tic;
format long
alpha=1.1;
nu=1;
k1=5000;
delta_r1(1)=0.05;
r(1)=0.3;
i=1;
while r(i)<0.5
    delta_r1(i+1)=alpha*delta_r1(i);
    r(i+1)=r(i)+delta_r1(i);
   
%     delta_r(i)=abs(r(i+1)-r(i));
        i=i+1;
end
Imax=i;
r(Imax)=0.5;
delta_r1(Imax-1)=0.5-r(Imax-1);
for i=1:Imax-1
delta_r(i)=delta_r1(i);
end


delta_teta=2*pi/Imax;
teta=[0:delta_teta:2*pi-delta_teta];
w=ones(Imax,Imax,k1);
w(:,1,:)=-5.625;
w(:,Imax,:)=-5.625;
 
for j=1:Imax-1
    A(j)=(2)/((delta_r(j)^2)*alpha*(alpha+1));
    B(j)=1/(r(j)*alpha*delta_r(j)*(1+alpha));
    C(j)=1/(r(j)*delta_teta)^2;
    %E(j)=1/(r(j)*A(j))+2/(A(j));
    %F(j)=-alpha^2*(1/(r(j)*A(j)))+2*alpha/(A(j)*delta_r(j));
end
  
D=zeros(Imax,Imax,k1);
si=ones(Imax,Imax,k1);
si(:,Imax,:)=1;
si1=ones(Imax,Imax,k1);
si0=ones(Imax,Imax,k1);

v=1;
k=1;
while v>10^(-6)
%    for k=1:1000   
      for i=1:Imax
si0(i,1,k+1)=si(i,2,k+1)-(2*(alpha)*delta_r(1)*1.5);
end

for i=2:Imax-1
        
        si(i,1,k+1)=(-si(i,2,k+1)*(A(1)+B(1))-si0(i,1,k)*(A(1)*alpha-(alpha^2)*B(1))-si(i+1,1,k+1)*C(1)-si(i-1,1,k)*C(1)-w(i,1,k))/(-A(1)*(alpha+1)+B(1)*((alpha^2)-1)-2*C(1));

     
end
        si(Imax,1,k+1)=(-si(Imax,2,k+1)*(A(1)+B(1))-si0(Imax,1,k)*(A(1)*alpha-(alpha^2)*B(1))-si(1,1,k+1)*C(1)-si(Imax-1,1,k)*C(1)-w(Imax,1,k))/(-A(1)*(alpha+1)+B(1)*((alpha^2)-1)-2*C(1));
        si(1,1,k+1)=(-si(1,2,k+1)*(A(1)+B(1))-si0(1,1,k)*(A(1)*alpha-(alpha^2)*B(1))-si(2,1,k+1)*C(1)-si(Imax,1,k)*C(1)-w(1,1,k))/(-A(1)*(alpha+1)+B(1)*((alpha^2)-1)-2*C(1));
% si(:,1,k)=(3*delta_r(1)-4*si(1,2,k)-si(1,3,k))/3;

for j=2:Imax-1
    
    for i=2:Imax-1
        
        si(i,j,k+1)=(-si(i,j+1,k+1)*(A(j)+B(j))-si(i,j-1,k)*(A(j)*alpha-(alpha^2)*B(j))-si(i+1,j,k+1)*C(j)-si(i-1,j,k)*C(j)-w(i,j,k))/(-A(j)*(alpha+1)+B(j)*((alpha^2)-1)-2*C(j));

     
    end
        si(Imax,j,k+1)=(-si(Imax,j+1,k+1)*(A(j)+B(j))-si(Imax,j-1,k)*(A(j)*alpha-(alpha^2)*B(j))-si(1,j,k+1)*C(j)-si(Imax-1,j,k)*C(j)-w(Imax,j,k))/(-A(j)*(alpha+1)+B(j)*((alpha^2)-1)-2*C(j));
 
        si(1,j,k+1)=(-si(1,j+1,k+1)*(A(j)+B(j))-si(1,j-1,k)*(A(j)*alpha-(alpha^2)*B(j))-si(2,j,k+1)*C(j)-si(Imax,j,k)*C(j)-w(1,j,k))/(-A(j)*(alpha+1)+B(j)*((alpha^2)-1)-2*C(j));
end

%  for i=2:Imax-1
%    si(i,1,k+1)=si(i,2,k+1)-r(1)*5*delta_r(1);
% 
%  end
%  si(Imax,1,k+1)=si(Imax,2,k+1)-r(1)*5*delta_r(1);
%  si(1,1,k+1)=si(1,2,k+1)-r(1)*5*delta_r(1);
  v=max(max(abs(si(:,:,k+1)-si(:,:,k))));
 
for j=2:Imax-1
    for i=2:Imax-1
        D(1,j,k+1)=nu*(alpha+1)*A(j)-((alpha^2)-1)*nu*B(j)+2*nu*C(j)+((alpha^2)-1)*B(j)*(si(2,j,k+1)-si(Imax,j,k+1))/(2*delta_teta);

    
        D(i,j,k+1)=nu*(alpha+1)*A(j)-((alpha^2)-1)*nu*B(j)+2*nu*C(j)+((alpha^2)-1)*B(j)*(si(i+1,j,k+1)-si(i-1,j,k+1))/(2*delta_teta);

    end
        D(Imax,j,k+1)=nu*(alpha+1)*A(j)-((alpha^2)-1)*nu*B(j)+2*nu*C(j)+((alpha^2)-1)*B(j)*(si(1,j,k+1)-si(Imax-1,j,k+1))/(2*delta_teta);

end

   for j=2:Imax-1

      for i=2:Imax-1
        w(1,j,k+1)=(w(1,j+1,k+1)*(nu*A(j)+nu*B(j)-B(j)*((si(2,j,k+1)-si(Imax,j,k+1))/(2*delta_teta)))+...
            w(1,j-1,k)*(alpha*nu*A(j)-(alpha^2)*nu*B(j)+(alpha^2)*B(j)*(si(2,j,k+1)-si(Imax-1,j,k+1))/(2*delta_teta))+...
            w(2,j,k+1)*(nu*C(j)-B(j)*(si(1,j+1,k+1)-(alpha^2)*si(1,j-1,k+1)+((alpha^2)-1)*si(1,j,k+1))/(2*delta_teta))+...
            w(Imax,j,k)*(nu*C(j)+B(j)*(si(1,j+1,k+1)-(alpha^2)*si(1,j-1,k+1)+((alpha^2)-1)*si(1,j,k+1))/(2*delta_teta)))/D(1,j,k+1);
        
        w(i,j,k+1)=(w(i,j+1,k+1)*(nu*A(j)+nu*B(j)-B(j)*((si(i+1,j,k+1)-si(i-1,j,k+1))/(2*delta_teta)))+...
            w(i,j-1,k)*(alpha*nu*A(j)-(alpha^2)*nu*B(j)+(alpha^2)*B(j)*(si(i+1,j,k+1)-si(i-1,j,k+1))/(2*delta_teta))+...
            w(i+1,j,k+1)*(nu*C(j)-B(j)*(si(i,j+1,k+1)-(alpha^2)*si(i,j-1,k+1)+((alpha^2)-1)*si(i,j,k+1))/(2*delta_teta))+...
            w(i-1,j,k)*(nu*C(j)+B(j)*(si(i,j+1,k+1)-(alpha^2)*si(i,j-1,k+1)+((alpha^2)-1)*si(i,j,k+1))/(2*delta_teta)))/D(i,j,k+1);
         
         w(Imax,j,k+1)=(w(Imax,j+1,k+1)*(nu*A(j)+nu*B(j)-B(j)*((si(1,j,k+1)-si(Imax-1,j,k+1))/(2*delta_teta)))+...
            w(Imax,j-1,k)*(alpha*nu*A(j)-(alpha^2)*nu*B(j)+(alpha^2)*B(j)*(si(1,j,k+1)-si(Imax-1,j,k+1))/(2*delta_teta))+...
            w(1,j,k+1)*(nu*C(j)-B(j)*(si(1,j+1,k+1)-(alpha^2)*si(1,j-1,k+1)+((alpha^2)-1)*si(1,j,k+1))/(2*delta_teta))+...
            w(Imax-1,j,k)*(nu*C(j)+B(j)*(si(Imax,j+1,k+1)-(alpha^2)*si(Imax,j-1,k+1)+((alpha^2)-1)*si(Imax,j,k+1))/(2*delta_teta)))/D(Imax,j,k+1);
      end

  end
 k=k+1;
  end
%  
si1=si(:,:,k);
w1=w(:,:,k);
% 
% % [T R]=pol2cart(teta,r(1,:));
% % [X Y]= meshgrid(T,R); 
% % 
% % figure(1)
% % contour(X,Y,si1,30)
% % colorbar;
% % 
% % figure(2)
% % contour(X,Y,w1)
% % colorbar
% 
% %%
figure(1)
[x,y]=meshgrid(r,teta);
plot(si1(1,:),r)

% figure(1)
[b,o]=pol2cart(y,x);
% contour(y,x,si1)
% 
figure(2)
contour(b,o,si1,30)

figure(3)
%  plot(w1(1,:),r)
 contour(b,o,w1,30)

% %%
% for j=1:Imax
%     for i=2:Imax-1
% Vr(i,j)=(si(i+1,j,k)-si(i-1,j,k))/(r(j)*2*delta_teta);
%     end
% Vr(1,j)=(si(2,j,k)-si(Imax,j,k))/(r(j)*2*delta_teta);
% Vr(Imax,j)=(si(1,j,k)-si(Imax-1,j,k))/(r(j)*2*delta_teta);
% 
% end
% for j=1:Imax-1
%     for i=1:Imax-1
% Vt(i,j)=(si(i+1,j,k+1)-si(i,j,k))/(2*delta_r(j));
%     end
% Vt(Imax,j)=(si(1,j,k+1)-si(Imax,j,k))/(2*delta_r(j));
% 
% end
% figure(8)
% plot(r,Vt);
toc;


