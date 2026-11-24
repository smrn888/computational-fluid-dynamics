clc
clear all
alpha=1;
nu=1;
k1=101;
wi=5;
delta_r(1)=0.1;
r(1)=0.3;
% if alpha<1
    i=1;
    while r(i)<0.5
    delta_r(i+1)=alpha*delta_r(i);
    r(i+1)=r(i)+delta_r(i);
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
delta_teta=ones(Imax,1)*2*pi/Imax;
teta=[0:delta_teta:2*pi-delta_teta];
w=zeros(Imax,Imax,1);
w(:,1,1)=-5.625;
w(:,Imax,1)=-5.625;

  for i=1:Imax
    for j=1:Imax
        A(j)=(alpha+alpha^2)*delta_r(j);
        B(i,j)=(alpha^2-1)/(r(j)*A(j))-2*(alpha+1)/(delta_r(j)*A(j))-2/((r(j)*delta_teta(i))^2);
        C(i,j)=1/(r(j)*delta_teta(i))^2;
        
    end
  end   
  D=zeros(Imax,Imax,k1);
  si=ones(Imax,Imax,k1);
 %si0=zeros(Imax,1,k1);
% for k=1:100
%     si(:,Imax,k)=0;

for i=2:Imax-1

    for j=2:Imax-1
        
%             si(1,1,1)=(((-si(1,2,1)*((1/(r(1)*A(1)))+2/(A(1)*delta_r(1))))+si0(1,1,1)*((alpha^2/(A(1)*r(1)))-(2*alpha/(A(1)*delta_r(1))))-si(2,1,1)*C(1,1)-si(Imax,1,1)*C(1,1) )-w(1,1))/B(1,1);

%      si(i,1,1)=(((-si(i,2,1)*(1/(r(1)*A(1))+2/(A(1)*delta_r(1))))+si0(i,1,1)*((alpha^2/(A(1)*r(1)))-(2*alpha/(A(1)*delta_r(1))))-si(i+1,1,1)*C(i,1)-si(i-1,1,1)*C(i,1) )-w(i,1))/B(i,1);

%         si(i,j,1)=((-si(i,j+1,1)*(1/(r(j)*A(j))+2/(A(j)*delta_r(j))))...
%             +si(i,j-1,1)*((alpha^2/(A(j)*r(j)))-(2*alpha/(A(j)*delta_r(j))))...
%             -si(i+1,j,1)*C(i,j)-si(i-1,j,1)*C(i,j)-w(i,j))/B(i,j);
        si(1,j,1)=(((-si(1,j+1,1)*((1/(r(j)*A(j)))+2/(A(j)*delta_r(j))))+si(1,j-1,1)*((alpha^2/(A(j)*r(j)))-(2*alpha/(A(j)*delta_r(j))))-si(2,j,1)*C(1,j)-si(Imax,j,1)*C(1,j) )-w(1,j))/B(1,j);
        si(i,j,1)=   (   -si(i,j+1,1)*(    1/(r(j)*A(j))  +   2/(A(j)*delta_r(j))     )+ si(i,j-1,1)*(    alpha^2/(A(j)*r(j))  -   2*alpha/(A(j)*delta_r(j))    ) -si(i+1,j,1)*C(1,j) -si(i-1,j,1)*C(1,j) -w(i,j)  )  /B(i,j);
                si(2,2,1)=   (   -si(2,3,1)*(    1/(r(2)*A(2))  +   2/(A(2)*delta_r(2))     )+ si(2,1,1)*(    alpha^2/(A(2)*r(2))  -   2*alpha/(A(2)*delta_r(2))    ) -si(3,j,1)*C(2,2) -si(1,2,1)*C(2,2) -w(2,2)  )  /B(2,2);

        si(Imax,j,1)=(((-si(Imax,j+1,1)*(1/(r(j)*A(j))+2/(A(j)*delta_r(j))))+si(Imax,j-1,1)*((alpha^2/(A(j)*r(j)))-(2*alpha/(A(j)*delta_r(j))))-si(1,j,1)*C(Imax,j)-si(Imax-1,j,1)*C(Imax,j) )-w(Imax,j))/B(Imax,j);
%         
        
        
%          si0(1,1,1)=si(2,1,1)-(1+1/alpha)*delta_r(1)*r(1)*wi;
%          si0(Imax,1,1)=si(1,1,1)-(1+1/alpha)*delta_r(1)*r(1)*wi;
%          si0(i,1,1)=si(i+1,1,1)-(1+1/alpha)*delta_r(1)*r(1)*wi;
%            si(i,1,1)=(((-si(i,2,1)*(1/(r(1)*A(1))+2/(A(1)*delta_r(1))))+si0(i,1,1)*((alpha^2/(A(1)*r(1)))-(2*alpha/(A(1)*delta_r(1))))-si(i+1,1,1)*C(i,1)-si(i-1,1,1)*C(i,1) )-w(i,1))/B(i,1);

%         si(i,j,k+1)=((-si(i,j+1,k)*(1/(r(j)*A(j)))+si(i,j-1,k+1)*(alpha^2/(A(j)*r(j)))-si(i+1,j,k)*C(i,j)-si(i-1,j,k+1)*C(i,j) )-w(i,j))/B(i,j);
%         si(1,j,k+1)=((-si(1,j+1,k)*(1/(r(j)*A(j)))+si(1,j-1,1)*(alpha^2/(A(j)*r(j)))-si(2,j,1)*C(1,j)-si(Imax,j,1)*C(1,j) )-w(1,j))/B(1,j);
%         si(Imax,j,1)=((-si(Imax,j+1,1)*(1/(r(j)*A(j)))+si(Imax,j-1,1)*(alpha^2/(A(j)*r(j)))-si(1,j,1)*C(Imax,j)-si(Imax-1,j,1)*C(Imax,j) )-w(Imax,j))/B(Imax,j);
%         D(i,j,1)=(si(i+1,j,1)-si(i-1,j,1))*(alpha-1)/(r(j)*2*delta_teta(j)*alpha*delta_r(j))+2*nu/(alpha*delta_r(j)^2)-nu*(alpha-1)/(r(j)*alpha*delta_r(j))+2*nu/(r(j)^2*delta_teta(j)^2);
%         D(1,j,1)=(si(2,j,1)-si(Imax,j,1))*(alpha-1)/(r(j)*2*delta_teta(j)*alpha*delta_r(j))+2*nu/(alpha*delta_r(j)^2)-nu*(alpha-1)/(r(j)*alpha*delta_r(j))+2*nu/(r(j)^2*delta_teta(j)^2);
%         D(Imax,j,1)=(si(1,j,1)-si(Imax-1,j,1))*(alpha-1)/(r(j)*2*delta_teta(j)*alpha*delta_r(j))+2*nu/(alpha*delta_r(j)^2)-nu*(alpha-1)/(r(j)*alpha*delta_r(j))+2*nu/(r(j)^2*delta_teta(j)^2);
%         w(i,j,1)=((si(i-1,j,1)-si(i+1,j,1))*(w(i,j+1,1)-alpha^2*w(i,j-1,1))/(r(j)*2*delta_teta(j)*A(j))+...
%             p((alpha^2-1)*si(i,j,1)+si(i,j+1,1)-alpha^2*si(i,j-1,1))*(w(i+1,j,1)-w(i-1,j,1))/(r(j)*2*delta_teta(j)*A(j))+...
%             nu*((w(i,j+1,1)+alpha*w(i,j-1,1))/(A(j)*delta_r(j)/2)+(w(i,j+1,1)-alpha^2*w(i,j-1,1))/(r(j)*A(j))+(w(i+1,j,1)+w(i-1,j,1))/C(j)))/(D(i,j,1));
%         w(1,j,1)=((si(Imax,j,1)-si(2,j,1))*(w(1,j+1,1)-alpha^2*w(1,j-1,1))/(r(j)*2*delta_teta(j)*A(j))+...
%             ((alpha^2-1)*si(1,j,1)+si(1,j+1,1)-alpha^2*si(1,j-1,1))*(w(2,j,1)-w(Imax,j,1))/(r(j)*2*delta_teta(j)*A(j))+...
%             nu*((w(1,j+1,1)+alpha*w(1,j-1,1))/(A(j)*delta_r(j)/2)+(w(1,j+1,1)-alpha^2*w(1,j-1,1))/(r(j)*A(j))+(w(2,j,1)+w(Imax,j,1))/C(j)))/(D(1,j,1));
%        w(Imax,j,1)=((si(Imax-1,j,1)-si(1,j,1))*(w(Imax,j+1,1)-alpha^2*w(Imax,j-1,1))/(r(j)*2*delta_teta(j)*A(j))+...
%             ((alpha^2-1)*si(Imax,j,1)+si(Imax,j+1,1)-alpha^2*si(Imax,j-1,1))*(w(1,j,1)-w(Imax-1,j,1))/(r(j)*2*delta_teta(j)*A(j))+...
%             nu*((w(Imax,j+1,1)+alpha*w(Imax,j-1,1))/(A(j)*delta_r(j)/2)+(w(Imax,j+1,1)-alpha^2*w(Imax,j-1,1))/(r(j)*A(j))+(w(1,j,1)+w(Imax-1,j,1))/C(j)))/(D(i,j,1));
%             
    end
%     si(Imax,1,1)=(((-si(Imax,2,1)*(1/(r(1)*A(1))+2/(A(1)*delta_r(1))))+si0(Imax,1,1)*((alpha^2/(A(1)*r(1)))-(2*alpha/(A(1)*delta_r(1))))-si(1,1,1)*C(Imax,1)-si(Imax-1,1,1)*C(Imax,1) )-w(Imax,1))/B(Imax,1);

%     si(i,1,1)=si(i,2,1)+r(1)*w(i,1)*delta_r(1);
%     si(1,1,1)=si(Imax,1,1)+r(1)*w(Imax,1)*delta_r(1);
%     si(Imax,1,1)=si(Imax,2,1)+r(1)*w(i,1)*delta_r(1);

end
% end

for i=2:Imax-1
               w(1,1)=-si(1,1,1)*B(1,1)+((-si(1,2,1)*(1/(r(1)*A(1))+2/(A(1)*delta_r(1))))+si0(1,1,1)*((alpha^2/(A(1)*r(1)))-(2*alpha/(A(1)*delta_r(1))))-si(2,1,1)*C(1,1)-si(Imax,1,1)*C(1,1));
    w(i,1)=-si(i,1,1)*B(i,1)+((-si(i,2,1)*(1/(r(1)*A(1))+2/(A(1)*delta_r(1))))+si0(i,1,1)*((alpha^2/(A(1)*r(1)))-(2*alpha/(A(1)*delta_r(1))))-si(i+1,1,1)*C(i,1)-si(i-1,1,1)*C(i,1));
        w(Imax,1)=-si(Imax,1,1)*B(Imax,1)+((-si(Imax,2,1)*(1/(r(1)*A(1))+2/(A(1)*delta_r(1))))+si0(Imax,1,1)*((alpha^2/(A(1)*r(1)))-(2*alpha/(A(1)*delta_r(1))))-si(1,1,1)*C(Imax,1)-si(Imax-1,1,1)*C(Imax,1));

end

