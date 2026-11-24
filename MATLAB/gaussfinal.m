clc
clear all
alpha=1.1;
nu=1;
k1=1000;
delta_r(1)=0.05;
r(1)=0.3;
% if alpha<1
    i=1;
    while r(i)<0.5
    delta_r(i+1)=alpha*delta_r(i);
    r(i+1)=r(i)+delta_r(i);
    if r(i)>0.5
        r(i)=0.5;
    end
%     delta_r(i)=abs(r(i+1)-r(i));
        i=i+1;
    end
       Imax=i;
%         r(Imax+1)=0.5;
%         delta_r(Imax)=0.5-r(Imax);
%  
% end
% if alpha>1
%     i=1;
%     while delta_r(i)>0.001
%     delta_r(i+1)=alpha*delta_r(i);
%     i=i+1;
%     end
% end

% 
% u=zeros
% for i=1:100
delta_teta=2*pi/Imax;
teta=[0:delta_teta:2*pi-delta_teta];
w=ones(Imax,Imax,k1);
 w(:,1,:)=-5.625;
 w(:,Imax,:)=-5.625;
 
    for j=1:Imax
        A(j)=(alpha+alpha^2)*delta_r(j);
        B(j)=(alpha^2-1)/(r(j)*A(j))-2*(alpha+1)/(delta_r(j)*A(j))-2/((r(j)*delta_teta)^2);
        C(j)=1/(r(j)*delta_teta)^2;
        E(j)=1/(r(j)*A(j))+2/(A(j));
        F(j)=-alpha^2*(1/(r(j)*A(j)))+2*alpha/(A(j)*delta_r(j));
      
        
    end
     
  D=zeros(Imax,Imax,k1);
  si=zeros(Imax,Imax,k1);
  si(:,Imax,:)=1;
  si1=zeros(Imax,Imax,k1);

  
v=1;
k=1;
while v>10^(-6)
%  for k=1:100    
for j=2:Imax-1
    
        si(1,j,k+1)=(-si(1,j+1,k)*E(j)-si(1,j-1,k+1)*F(j)-si(2,j,k)*C(j)-si(Imax,j,k+1)*C(j)-w(1,j,k))/B(j);

    for i=2:Imax-1
        
        si(i,j,k+1)=(-si(i,j+1,k)*E(j)-si(i,j-1,k+1)*F(j)-si(i+1,j,k)*C(j)-si(i-1,j,k+1)*C(j)-w(i,j,k))/B(j);

     
    end
        si(Imax,j,k+1)=(-si(Imax,j+1,k)*E(j)-si(Imax,j-1,k+1)*F(j)-si(1,j,k)*C(j)-si(Imax-1,j,k+1)*C(j)-w(Imax,j,k))/B(j);
 

end
 for i=2:Imax-1
   si(i,1,k+1)=si(i,2,k+1)+r(1)*5*delta_r(1);

 end
   si(Imax,1,k+1)=si(Imax,2,k+1)+r(1)*5*delta_r(1);
   si(1,1,k+1)=si(1,2,k+1)+r(1)*5*delta_r(1);
 v=max(max(abs(si(:,:,k+1)-si(:,:,k))));
 
for j=2:Imax-1
    for i=2:Imax-1
        D(1,j,k+1)=(((alpha^2)-1)*(si(2,j,k+1)-si(Imax,j,k+1))/(A(j)*2*r(j)*delta_teta))-nu*B(j);

    
        D(i,j,k+1)=(((alpha^2)-1)*(si(i+1,j,k+1)-si(i-1,j,k+1))/(A(j)*2*r(j)*delta_teta))-nu*B(j);

    end
        D(Imax,j,k+1)=(((alpha^2)-1)*(si(1,j,k+1)-si(Imax-1,j,k+1))/(A(j)*2*r(j)*delta_teta))-nu*B(j);

end

   for i=2:Imax-1

      for j=2:Imax-1
        w(1,j,k+1)=(w(1,j+1,k)*(((nu/(A(j)*r(j)))+2*nu/(A(j)*delta_r(j)))-(si(2,j,k+1)-si(Imax,j,k+1))/(2*r(j)*A(j)*delta_teta))+...
            w(1,j-1,k+1)*(-(alpha^2)*nu/(A(j)*r(j))+2*nu*alpha/(delta_r(j)*A(j))+(alpha^2)*(si(2,j,k+1)-si(Imax,j,k+1))/(2*r(j)*delta_teta*A(j)))+...
            w(2,j,k)*(nu*C(j)+(si(1,j+1,k+1)-(alpha^2)*si(1,j-1,k+1)+((alpha^2)-1)*si(1,j,k+1))/(2*r(j)*delta_teta*A(j)))+...
            w(Imax,j,k+1)*(nu*C(j)-(si(1,j+1,k+1)-(alpha^2)*si(1,j-1,k+1)+((alpha^2)-1)*si(1,j,k+1))/(2*r(j)*delta_teta*A(j))))/D(1,j,k+1);
        
        w(i,j,k+1)=(w(i,j+1,k)*(((nu/(A(j)*r(j)))+2*nu/(A(j)*delta_r(j)))-(si(i+1,j,k+1)-si(i-1,j,k+1))/(2*r(j)*A(j)*delta_teta))+...
            w(i,j-1,k+1)*(-(alpha^2)*nu/(A(j)*r(j))+2*nu*alpha/(delta_r(j)*A(j))+(alpha^2)*(si(i+1,j,k+1)-si(i-1,j,k+1))/(2*r(j)*delta_teta*A(j)))+...
            w(i+1,j,k)*(nu*C(j)+(si(i,j+1,k+1)-(alpha^2)*si(i,j-1,k+1)+((alpha^2)-1)*si(i,j,k+1))/(2*r(j)*delta_teta*A(j)))+...
            w(i-1,j,k+1)*(nu*C(j)-(si(i,j+1,k+1)-(alpha^2)*si(i,j-1,k+1)+((alpha^2)-1)*si(i,j,k+1))/(2*r(j)*delta_teta*A(j))))/D(i,j,k+1);
         
         w(Imax,j,k+1)=(w(Imax,j+1,k)*(((nu/(A(j)*r(j)))+2*nu/(A(j)*delta_r(j)))-(si(1,j,k+1)-si(Imax-1,j,k+1))/(2*r(j)*A(j)*delta_teta))+...
            w(Imax,j-1,k+1)*(-(alpha^2)*nu/(A(j)*r(j))+2*nu*alpha/(delta_r(j)*A(j))+(alpha^2)*(si(1,j,k+1)-si(Imax-1,j,k+1))/(2*r(j)*delta_teta*A(j)))+...
            w(1,j,k)*(nu*C(j)+(si(Imax,j+1,k+1)-(alpha^2)*si(Imax,j-1,k+1)+((alpha^2)-1)*si(Imax,j,k+1))/(2*r(j)*delta_teta*A(j)))+...
            w(Imax-1,j,k+1)*(nu*C(j)-(si(Imax,j+1,k+1)-(alpha^2)*si(Imax,j-1,k+1)+((alpha^2)-1)*si(Imax,j,k+1))/(2*r(j)*delta_teta*A(j))))/D(Imax,j,k+1);    
      end

  end
 k=k+1;
end
si1=si(:,:,k);
w1=w(:,:,k);
% 
% [T R]=pol2cart(teta,r(1,:));
% [X Y]= meshgrid(T,R); 
% 
% figure(1)
% contour(X,Y,si1,30)
% colorbar;
% 
% figure(2)
% contour(X,Y,w1)
% colorbar

%%
[x,y]=meshgrid(r,teta);
[b,o]=pol2cart(y,x);
contour(b,o,si1,100)
%%
figure
contour(b,o,w1,100)
%%
% for j=1:Imax
%     for i=2:Imax-1
% Vr(i,j)=(si(i+1,j,k)-si(i-1,j))/(r(j)*2*delta_teta);
%     end
% Vr(1,j)=(si(2,j,k)-si(Imax,j))/(r(j)*2*delta_teta);
% Vr(Imax,j)=(si(1,j,k)-si(Imax-1,j))/(r(j)*2*delta_teta);
% 
% end
% for j=1:Imax
%     for i=1:Imax-1
% Vt(i,j)=(si(i+1,j,k)-si(i,j))/(2*delta_r(j));
%     end
% Vt(Imax,j)=(si(1,j,k)-si(Imax,j))/(2*delta_r(j));
% 
% end
% figure
% plot(Vt,r);
% 

