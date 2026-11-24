clc
clear all
alpha=1;
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
delta_teta=ones(Imax,1)*2*pi/Imax;
teta=[0:delta_teta:2*pi-delta_teta];
w=ones(Imax,Imax,k1);
 w(:,1,:)=-5.625;
 w(:,Imax,:)=-5.625;
 
  for i=1:Imax
    for j=1:Imax
        A(j)=(alpha+alpha^2)*delta_r(j);
        B(i,j)=(alpha^2-1)/(r(j)*A(j))-2*(alpha+1)/(delta_r(j)*A(j))-2/((r(j)*delta_teta(i))^2);
        C(i,j)=1/(r(j)*delta_teta(i))^2;
        E(i,j)=1/(r(j)*A(j))+2/(A(j));
        F(i,j)=-alpha^2*(1/(r(j)*A(j)))+2*alpha/(A(j)*delta_r(j));
      
        
    end
  end   
  D=zeros(Imax,Imax,k1);
  si=ones(Imax,Imax,k1);
  si1=zeros(Imax,Imax,k1);

  
v=1;
k=1;
while v>10^(-12)
%  for k=1:100    
for j=2:Imax-1

    for i=2:Imax-1
        
        si(i,j,k+1)=(-si(i,j+1,k)*E(i,j)-si(i,j-1,k+1)*F(i,j)-si(i+1,j,k)*C(i,j)-si(i-1,j,k+1)*C(i,j)-w(i,j,k))/B(i,j);

     
    end
        si(1,j,k+1)=(-si(1,j+1,k)*E(1,j)-si(1,j-1,k+1)*F(1,j)-si(2,j,k+1)*C(1,j)-si(Imax,j,k+1)*C(1,j)-w(1,j,k))/B(1,j);
        si(Imax,j,k+1)=(-si(Imax,j+1,k)*E(Imax,j)-si(Imax,j-1,k+1)*F(Imax,j)-si(1,j,k+1)*C(Imax,j)-si(Imax-1,j,k+1)*C(Imax,j)-w(Imax,j,k))/B(Imax,j);
 

end
 for i=2:Imax-1
   si(i,1,k+1)=si(i,2,k+1)+r(1)*5*delta_r(1);

 end
   si(Imax,1,k+1)=si(Imax,2,k+1)+r(1)*5*delta_r(1);
   si(1,1,k+1)=si(1,2,k+1)+r(1)*5*delta_r(1);
 v=max(max(abs(si(:,:,k+1)-si(:,:,k))));

for i=2:Imax-1

    for j=2:Imax-1
        D(i,j,k+1)=(((alpha^2)-1)*(si(i+1,j,k+1)-si(i-1,j,k+1))/(A(j)*2*r(j)*delta_teta(j)))-nu*B(i,j);
        D(1,j,k+1)=(((alpha^2)-1)*(si(2,j,k+1)-si(Imax,j,k+1))/(A(j)*2*r(j)*delta_teta(j)))-nu*B(1,j);
        D(i,1,k+1)=(((alpha^2)-1)*(si(i+1,1,k+1)-si(i-1,1,k+1))/(A(1)*2*r(1)*delta_teta(1)))-nu*B(i,1);
        D(Imax,j,k+1)=(((alpha^2)-1)*(si(1,j,k+1)-si(Imax-1,j,k+1))/(A(j)*2*r(j)*delta_teta(j)))-nu*B(Imax,j);
        D(i,Imax,k+1)=(((alpha^2)-1)*(si(i+1,Imax,k+1)-si(i-1,Imax,k+1))/(A(Imax)*2*r(Imax)*delta_teta(Imax)))-nu*B(i,Imax);
        D(1,Imax,k+1)=(((alpha^2)-1)*(si(2,Imax,k+1)-si(Imax,Imax,k+1))/(A(Imax)*2*r(Imax)*delta_teta(Imax)))-nu*B(1,Imax);
        D(Imax,1,k+1)=(((alpha^2)-1)*(si(1,1,k+1)-si(Imax-1,1,k+1))/(A(1)*2*r(1)*delta_teta(1)))-nu*B(Imax,1);
        D(1,1,k+1)=(((alpha^2)-1)*(si(2,1,k+1)-si(Imax,1,k+1))/(A(1)*2*r(1)*delta_teta(1)))-nu*B(1,1);
        D(Imax,Imax,k+1)=(((alpha^2)-1)*(si(1,Imax,k+1)-si(Imax-1,Imax,k+1))/(A(Imax)*2*r(Imax)*delta_teta(Imax)))-nu*B(Imax,Imax);

    end
end

for i=2:Imax-1

    for j=2:Imax-1
       w(1,j,k+1)=(w(1,j+1,k)*(((nu/(A(j)*r(j)))+2*nu/(A(j)*delta_r(j)))-(si(2,j,k+1)-si(Imax,j,k+1))/(2*r(j)*A(j)*delta_teta(j)))+...
            w(1,j-1,k)*(-(alpha^2)*nu/(A(j)*r(j))+2*nu*alpha/(delta_r(j)*A(j))+(alpha^2)*(si(2,j,k+1)-si(Imax,j,k+1))/(2*r(j)*delta_teta(j)*A(j))+...
            w(2,j,k)*(nu*C(1,j)+(si(1,j+1,k+1)-(alpha^2)*si(1,j-1,k+1)+((alpha^2)-1)*si(1,j,k+1)/(2*r(j)*delta_teta(j)*A(j))))+...
            w(Imax,j,k)*(nu*C(1,j)-(si(1,j+1,k+1)-(alpha^2)*si(1,j-1,k+1)+((alpha^2)-1)*si(1,j,k+1)/(2*r(j)*delta_teta(j)*A(j))))))/D(1,j,k+1);
        
        w(i,j,k+1)=(w(i,j+1,k)*(((nu/(A(j)*r(j)))+2*nu/(A(j)*delta_r(j)))-(si(i+1,j,k+1)-si(i-1,j,k+1))/(2*r(j)*A(j)*delta_teta(j)))+...
            w(i,j-1,k+1)*(-(alpha^2)*nu/(A(j)*r(j))+2*nu*alpha/(delta_r(j)*A(j))+(alpha^2)*(si(i+1,j,k+1)-si(i-1,j,k+1))/(2*r(j)*delta_teta(j)*A(j))+...
            w(i+1,j,k)*(nu*C(i,j)+(si(i,j+1,k+1)-(alpha^2)*si(i,j-1,k+1)+((alpha^2)-1)*si(i,j,k+1)/(2*r(j)*delta_teta(j)*A(j))))+...
            w(i-1,j,k)*(nu*C(i,j)-(si(i,j+1,k+1)-(alpha^2)*si(i,j-1,k+1)+((alpha^2)-1)*si(i,j,k+1)/(2*r(j)*delta_teta(j)*A(j))))))/D(i,j,k+1);
         
         w(Imax,j,k+1)=(w(Imax,j+1,k)*(((nu/(A(j)*r(j)))+2*nu/(A(j)*delta_r(j)))-(si(1,j,k+1)-si(Imax-1,j,k+1))/(2*r(j)*A(j)*delta_teta(j)))+...
            w(Imax,j-1)*(-(alpha^2)*nu/(A(j)*r(j))+2*nu*alpha/(delta_r(j)*A(j))+(alpha^2)*(si(1,j,k+1)-si(Imax-1,j,k+1))/(2*r(j)*delta_teta(j)*A(j))+...
            w(1,j,k)*(nu*C(Imax,j)+(si(Imax,j+1,k+1)-(alpha^2)*si(Imax,j-1,k+1)+((alpha^2)-1)*si(Imax,j,k+1)/(2*r(j)*delta_teta(j)*A(j))))+...
            w(Imax-1,j,k)*(nu*C(Imax,j)-(si(Imax,j+1,k+1)-(alpha^2)*si(Imax,j-1,k+1)+((alpha^2)-1)*si(Imax,j,k+1)/(2*r(j)*delta_teta(j)*A(j))))))/D(Imax,j,k+1);
    end
end


 k=k+1;
end