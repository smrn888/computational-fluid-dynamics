clc
clear all
alpha=1;
nu=0.01;
k1=1000;
delta_r(1)=0.05;
r(1)=0.3;
    i=1;
    while r(i)<0.5
        
    delta_r(i+1)=alpha*delta_r(i);
    r(i+1)=r(i)+delta_r(i);
          i=i+1;
    end
      if r(i)>0.5
        r(i)=0.5;
      end
       Imax=i;
delta_r=delta_r(1:Imax-1);
delta_r(i-1)=r(i)-r(i-1);
delta_teta=2*pi/Imax;



w=ones(Imax,Imax);
 w(:,1)=-5.625;
 w(:,Imax)=-5.625;
 
    for j=1:Imax-1
        A(j)=(alpha+alpha^2)*delta_r(j);
        B(j+1)=(alpha^2-1)/(r(j+1)*A(j))-2*(alpha+1)/(delta_r(j)*A(j))-2/((r(j+1)*delta_teta)^2);
        C(j+1)=1/(r(j+1)*delta_teta)^2;
        E(j+1)=1/(r(j+1)*A(j))+2/(A(j));
        F(j+1)=-alpha^2*(1/(r(j+1)*A(j)))+2*alpha/(A(j)*delta_r(j));
      
        
    end
  D=zeros(Imax,Imax);
  si=ones(Imax,Imax);
  si1=zeros(Imax,Imax);

  
v=1;
k=1
while v>10^(-12)
    
si1=si;
for j=2:Imax-1

    for i=2:Imax-1
        
        si(i,j)=(-si(i,j+1)*E(j)-si(i,j-1)*F(j)-si(i+1,j)*C(j)-si(i-1,j)*C(j)-w(i,j))/B(j);

     
    end
        si(1,j)=(-si(1,j+1)*E(j)-si(1,j-1)*F(j)-si(2,j)*C(j)-si(Imax,j)*C(j)-w(1,j))/B(j);
        si(Imax,j)=(-si(Imax,j+1)*E(j)-si(Imax,j-1)*F(j)-si(1,j)*C(j)-si(Imax-1,j)*C(j)-w(Imax,j))/B(j);
 

end
 for i=2:Imax-1
   si(i,1)=si(i,2)+r(1)*5*delta_r(1);

 end
   si(Imax,1)=si(Imax,2)+r(1)*5*delta_r(1);
   si(1,1)=si(1,2)+r(1)*5*delta_r(1);
v=max(max(abs(si1-si)));

for i=2:Imax-1

    for j=2:Imax-1
        D(i,j)=(((alpha^2)-1)*(si(i+1,j)-si(i-1,j))/(A(j)*2*r(j)*delta_teta))-nu*B(j);
        D(1,j)=(((alpha^2)-1)*(si(2,j)-si(Imax,j))/(A(j)*2*r(j)*delta_teta))-nu*B(j);
        D(Imax,j)=(((alpha^2)-1)*(si(1,j)-si(Imax-1,j))/(A(j)*2*r(j)*delta_teta))-nu*B(j);

    end
end

for i=2:Imax-1

    for j=2:Imax-1
       w(1,j)=(w(1,j+1)*(((nu/(A(j)*r(j)))+2*nu/(A(j)*delta_r(j)))-(si(2,j)-si(Imax,j))/(2*r(j)*A(j)*delta_teta))+...
            w(1,j-1)*(-(alpha^2)*nu/(A(j)*r(j))+2*nu*alpha/(delta_r(j)*A(j))+(alpha^2)*(si(2,j)-si(Imax,j))/(2*r(j)*delta_teta*A(j))+...
            w(2,j)*(nu*C(j)+(si(1,j+1)-(alpha^2)*si(1,j-1)+((alpha^2)-1)*si(1,j)/(2*r(j)*delta_teta*A(j))))+...
            w(Imax,j)*(nu*C(j)-(si(1,j+1)-(alpha^2)*si(1,j-1)+((alpha^2)-1)*si(1,j)/(2*r(j)*delta_teta*A(j))))))/D(1,j);
        
        w(i,j)=(w(i,j+1)*(((nu/(A(j)*r(j)))+2*nu/(A(j)*delta_r(j)))-(si(i+1,j)-si(i-1,j))/(2*r(j)*A(j)*delta_teta))+...
            w(i,j-1)*(-(alpha^2)*nu/(A(j)*r(j))+2*nu*alpha/(delta_r(j)*A(j))+(alpha^2)*(si(i+1,j)-si(i-1,j))/(2*r(j)*delta_teta*A(j))+...
            w(i+1,j)*(nu*C(j)+(si(i,j+1)-(alpha^2)*si(i,j-1)+((alpha^2)-1)*si(i,j)/(2*r(j)*delta_teta*A(j))))+...
            w(i-1,j)*(nu*C(j)-(si(i,j+1)-(alpha^2)*si(i,j-1)+((alpha^2)-1)*si(i,j)/(2*r(j)*delta_teta*A(j))))))/D(i,j);
         
         w(Imax,j)=(w(Imax,j+1)*(((nu/(A(j)*r(j)))+2*nu/(A(j)*delta_r(j)))-(si(1,j)-si(Imax-1,j))/(2*r(j)*A(j)*delta_teta))+...
            w(Imax,j-1)*(-(alpha^2)*nu/(A(j)*r(j))+2*nu*alpha/(delta_r(j)*A(j))+(alpha^2)*(si(1,j)-si(Imax-1,j))/(2*r(j)*delta_teta*A(j))+...
            w(1,j)*(nu*C(Imax)+(si(Imax,j+1)-(alpha^2)*si(Imax,j-1)+((alpha^2)-1)*si(Imax,j)/(2*r(j)*delta_teta*A(j))))+...
            w(Imax-1,j)*(nu*C(Imax)-(si(Imax,j+1)-(alpha^2)*si(Imax,j-1)+((alpha^2)-1)*si(Imax,j)/(2*r(j)*delta_teta*A(j))))))/D(Imax,j);
    end
end


 k=k+1;
end