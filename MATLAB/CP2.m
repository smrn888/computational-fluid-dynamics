%% Getting Variable (delta_x) and Define mesh
clear all;
clc;
format long;

        delta_theta=input('Please enter delta theta by radian: ');
        delta_r =input('Please enter delta r for a course mesh: ');
        alfa=input('Please enter alfa: (delta x / alfa) is x steps: ');
        
        r=0.3;
        j=0;
        
        while (r<=0.5)
            
            r=r+delta_r/alfa^j;
            j=j+1;
            
        end
        
        j_max=j;
        matrix_dr =zeros(1,j_max-1);
        matrix_r=zeros(1,j_max);
        
        for k=1:j_max
            
            matrix_dr(1,k)=delta_r/alfa^k;
        
        end
        
          matrix_r(1,1)=0.3;
          matrix_r(1,j_max)=0.5;
       
          for k=1:j_max-1
            
              matrix_r(1,k+1)= matrix_r(1,k)+matrix_dr(1,k);
          end
        
          matrix_r(1,1)=0.3;
          matrix_r(1,j_max)=0.5;
        
        % Theta
        
          i_max = 2*pi/delta_theta; 
          matrix_dtheta =zeros(1,i_max);
          matrix_theta=zeros(1,i_max);
          matrix_theta(1,1)=0;
          matrix_theta(1,i_max)=2*pi;
        
          for k=1:i_max
              
             matrix_dtheta(1,k)=delta_theta;     
          
          end 
        
          for k=1:i_max-1
          
              matrix_theta(1,k+1)= matrix_theta(1,k)+ matrix_dtheta(1,k);
          
          end
          
t=input('for GS enter 1, for SOR enter 2, for SLOR enter 3, for multigrid enter 4, for SLOR+ADI enter 5:');
switch t
          
%% GAUSS-SEIDEL Method:
    case 1
k_max=1000;
        
fai=zeros(i_max,j_max,k_max);
sai=zeros(i_max,j_max,k_max);
e=zeros(i_max,j_max,k_max);
e(:,:,1)=1;

k=1;

%  while(e(:,:,k)>1e-10)
 for k=1:k_max
%SAI
sai(:,1,k)=(3*matrix_dr(1,1)-4*sai(1,2,k)-sai(1,3,k))/3;
for i=1:i_max
    
    if i==1
        
         A=2/(alfa*(alfa+1)*matrix_dr(1,1)^2)+(1/(matrix_r(1,2)*alfa*(alfa+1)*matrix_dr(1,1)));
         B=1/((matrix_r(1,2)^2)*delta_theta^2);
         C=2/((alfa+1)*matrix_dr(1,1)^2)-(alfa/(matrix_r(1,2)*(alfa+1)*matrix_dr(1,1)));
         D=-2/(alfa*matrix_dr(1,1)^2)-(2*B)+((alfa-1)/(matrix_r(1,2)*alfa*matrix_dr(1,1)));
         
         sai(1,2,k+1)=(-fai(1,2,k)+sai(1,3,k)*((C/3)-A)-sai(2,2,k)*B-(matrix_dr(1,1)*C)-sai(i_max,2,k)*B)/(D-(4*C/3));              
 
         
         for j=3:j_max-1
             
             
         A=2/(alfa*(alfa+1)*matrix_dr(1,j-1)^2)+(1/(matrix_r(1,j)*alfa*(alfa+1)*matrix_dr(1,j-1)));
         B=1/((matrix_r(1,j)^2)*delta_theta^2);
         C=2/((alfa+1)*matrix_dr(1,j-1)^2)-(alfa/(matrix_r(1,j)*(alfa+1)*matrix_dr(1,j-1)));
         D=-2/(alfa*matrix_dr(1,j-1)^2)-(2*B)+((alfa-1)/(matrix_r(1,j)*alfa*matrix_dr(1,j-1)));
         
         sai(1,j,k+1)=(-fai(1,j,k)+sai(1,j+1,k)*(-A)-sai(2,j,k)*B-sai(1,j-1,k+1)*C-sai(i_max,j,k)*B)/D;              
         
         end
    end
    
    if (i~=i_max) && (i~=1)
        
    for j=2:j_max-1
       
       if j==2
                
         A=2/(alfa*(alfa+1)*matrix_dr(1,1)^2)+(1/(matrix_r(1,2)*alfa*(alfa+1)*matrix_dr(1,1)));
         B=1/((matrix_r(1,2)^2)*delta_theta^2);
         C=2/((alfa+1)*matrix_dr(1,1)^2)-(alfa/(matrix_r(1,2)*(alfa+1)*matrix_dr(1,1)));
         D=-2/(alfa*matrix_dr(1,1)^2)-(2*B)+((alfa-1)/(matrix_r(1,2)*alfa*matrix_dr(1,1)));
         
         sai(i,2,k+1)=(-fai(i,2,k)+sai(i,3,k)*((C/3)-A)-sai(i+1,2,k)*B-(matrix_dr(1,1)*C)-sai(i-1,2,k+1)*B)/(D-(4*C/3));              
        
       end
       
           
         A=2/(alfa*(alfa+1)*matrix_dr(1,j-1)^2)+(1/(matrix_r(1,j)*alfa*(alfa+1)*matrix_dr(1,j-1)));
         B=1/((matrix_r(1,j)^2)*delta_theta^2);
         C=2/((alfa+1)*matrix_dr(1,j-1)^2)-(alfa/(matrix_r(1,j)*(alfa+1)*matrix_dr(1,j-1)));
         D=-2/(alfa*matrix_dr(1,j-1)^2)-(2*B)+((alfa-1)/(matrix_r(1,j)*alfa*matrix_dr(1,j-1)));
         
         sai(i,j,k+1)=(-fai(i,j,k)+sai(i,j+1,k)*(-A)-sai(i+1,j,k)*B-sai(i,j-1,k+1)*C-sai(i-1,j,k+1)*B)/D;              
         
    end
       
    end
    
    if i==i_max
       
         for j=2:j_max-1
       
       if j==2
                
         A=2/(alfa*(alfa+1)*matrix_dr(1,1)^2)+(1/(matrix_r(1,2)*alfa*(alfa+1)*matrix_dr(1,1)));
         B=1/((matrix_r(1,2)^2)*delta_theta^2);
         C=2/((alfa+1)*matrix_dr(1,1)^2)-(alfa/(matrix_r(1,2)*(alfa+1)*matrix_dr(1,1)));
         D=-2/(alfa*matrix_dr(1,1)^2)-(2*B)+((alfa-1)/(matrix_r(1,2)*alfa*matrix_dr(1,1)));
         
         sai(i_max,2,k+1)=(-fai(i_max,2,k)+sai(i_max,3,k)*((C/3)-A)-sai(1,2,k+1)*B-(matrix_dr(1,1)*C)-sai(i_max-1,2,k+1)*B)/(D-(4*C/3));              
        
       end
       
           
         A=2/(alfa*(alfa+1)*matrix_dr(1,j-1)^2)+(1/(matrix_r(1,j)*alfa*(alfa+1)*matrix_dr(1,j-1)));
         B=1/((matrix_r(1,j)^2)*delta_theta^2);
         C=2/((alfa+1)*matrix_dr(1,j-1)^2)-(alfa/(matrix_r(1,j)*(alfa+1)*matrix_dr(1,j-1)));
         D=-2/(alfa*matrix_dr(1,j-1)^2)-(2*B)+((alfa-1)/(matrix_r(1,j)*alfa*matrix_dr(1,j-1)));
         
         sai(i_max,j,k+1)=(-fai(i_max,j,k)+sai(i_max,j+1,k)*(-A)-sai(1,j,k)*B-sai(i_max,j-1,k+1)*C-sai(i_max-1,j,k+1)*B)/D;              
         
    
         end
  
    end
    
end

%FAI
fai(:,j_max,:)=-5.625;
fai(:,1,:)=-5.625;

for i=1:i_max
    
    mio=.0001;
    
    if i==1 
       
       P=(A-(1/mio)*(sai(i+1,2,k+1)-sai(i_max,2,k+1))/2*matrix_r(1,2)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,1));
       Q=(B+(1/mio)*(sai(i,3,k+1)+sai(i,2,k+1)*(alfa^2-1)-sai(i,1,k+1)*alfa^2)/2*matrix_r(1,2)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,2));
       R=((2/(alfa+1)*matrix_dr(1,1)^2)-alfa/(matrix_r(1,2)*(alfa+1)*matrix_dr(1,1))-(1/mio)*(-alfa*(sai(i+1,2,k+1)-sai(i_max,2,k+1)))/2*matrix_r(1,2)*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,1));
       S=B-(1/mio)*(sai(i,3,k+1)+(alfa^2-1)*sai(i,2,k+1)-(alfa^2)*sai(i,1,k+1))/(2*matrix_r(1,2)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,2));
       T=D-(1/mio)*(alfa-1)*(sai(i+1,2,k+1)-sai(i_max,2,k+1))/2*matrix_r(1,2)*(alfa)*matrix_dtheta(1,i)*matrix_dr(1,1);
       
       fai(i,2,k+1)=(-fai(i,3,k)*P-fai(i+1,2,k)*Q+5.625*R-fai(i_max,2,k)*S)/T;
       
       
         for j=3:j_max-1
             
       P=(A-(1/mio)*(sai(i+1,j,k+1)-sai(i_max,j,k+1))/2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j-1));
       Q=(B+(1/mio)*(sai(i,j+1,k+1)+sai(i,j,k+1)*(alfa^2-1)-sai(i,j-1,k+1)*alfa^2)/2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j));
       R=((2/(alfa+1)*matrix_dr(1,j-1)^2)-alfa/(matrix_r(1,j)*(alfa+1)*matrix_dr(1,j-1))-(1/mio)*(-alfa*(sai(i+1,j,k+1)-sai(i_max,j,k+1)))/2*matrix_r(1,j)*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j-1));
       S=B-(1/mio)*(sai(i,j+1,k+1)+(alfa^2-1)*sai(i,j,k+1)-(alfa^2)*sai(i,j-1,k+1))/(2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j));
       T=D-(1/mio)*(alfa-1)*(sai(i+1,j,k+1)-sai(i_max,j,k+1))/2*matrix_r(1,j)*(alfa)*matrix_dtheta(1,i)*matrix_dr(1,j-1);
       
       fai(i,j,k+1)=(-fai(i,j+1,k)*P-fai(i+1,j,k)*Q-fai(i,j-1,k+1)*R-fai(i_max,j,k)*S)/T; 
         
         end
    end
    
    if (i~=i_max) && (i~=1)
        
    for j=2:j_max-1
       
       if j==2
                      
       P=(A-(1/mio)*(sai(i+1,j,k+1)-sai(i-1,j,k+1))/2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j-1));
       Q=(B+(1/mio)*(sai(i,j+1,k+1)+sai(i,j,k+1)*(alfa^2-1)-sai(i,j-1,k+1)*alfa^2)/2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j));
       R=((2/(alfa+1)*matrix_dr(1,j-1)^2)-alfa/(matrix_r(1,j)*(alfa+1)*matrix_dr(1,j-1))-(1/mio)*(-alfa*(sai(i+1,j,k+1)-sai(i-1,j,k+1)))/2*matrix_r(1,j)*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j-1));
       S=B-(1/mio)*(sai(i,j+1,k+1)+(alfa^2-1)*sai(i,j,k+1)-(alfa^2)*sai(i,j-1,k+1))/(2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j));
       T=D-(1/mio)*(alfa-1)*(sai(i+1,j,k+1)-sai(i-1,j,k+1))/2*matrix_r(1,j)*(alfa)*matrix_dtheta(1,i)*matrix_dr(1,j-1);
       
       fai(i,j,k+1)=(-fai(i,j+1,k)*P-fai(i+1,j,k)*Q +5.625*R-fai(i-1,j,k+1)*S)/T; 
         
       end
       
       P=(A-(1/mio)*(sai(i+1,j,k+1)-sai(i-1,j,k+1))/2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j-1));
       Q=(B+(1/mio)*(sai(i,j+1,k+1)+sai(i,j,k+1)*(alfa^2-1)-sai(i,j-1,k+1)*alfa^2)/2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j));
       R=((2/(alfa+1)*matrix_dr(1,j-1)^2)-alfa/(matrix_r(1,j)*(alfa+1)*matrix_dr(1,j-1))-(1/mio)*(-alfa*(sai(i+1,j,k+1)-sai(i-1,j,k+1)))/2*matrix_r(1,j)*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j-1));
       S=B-(1/mio)*(sai(i,j+1,k+1)+(alfa^2-1)*sai(i,j,k+1)-(alfa^2)*sai(i,j-1,k+1))/(2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j));
       T=D-(1/mio)*(alfa-1)*(sai(i+1,j,k+1)-sai(i-1,j,k+1))/2*matrix_r(1,j)*(alfa)*matrix_dtheta(1,i)*matrix_dr(1,j-1);
       
       fai(i,j,k+1)=(-fai(i,j+1,k)*P-fai(i+1,j,k)*Q-fai(i,j-1,k+1)*R-fai(i-1,j,k+1)*S)/T; 
          
    end
       
    end
    
    if i==i_max
       
         for j=2:j_max-1
       
       if j==2
                
                      
       P=(A-(1/mio)*(sai(1,j,k+1)-sai(i-1,j,k+1))/2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j-1));
       Q=(B+(1/mio)*(sai(i,j+1,k+1)+sai(i,j,k+1)*(alfa^2-1)-sai(i,j-1,k+1)*alfa^2)/2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j));
       R=((2/(alfa+1)*matrix_dr(1,j-1)^2)-alfa/(matrix_r(1,j)*(alfa+1)*matrix_dr(1,j-1))-(1/mio)*(-alfa*(sai(1,j,k+1)-sai(i-1,j,k+1)))/2*matrix_r(1,j)*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j-1));
       S=B-(1/mio)*(sai(i,j+1,k+1)+(alfa^2-1)*sai(i,j,k+1)-(alfa^2)*sai(i,j-1,k+1))/(2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j));
       T=D-(1/mio)*(alfa-1)*(sai(1,j,k+1)-sai(i-1,j,k+1))/2*matrix_r(1,j)*(alfa)*matrix_dtheta(1,i)*matrix_dr(1,j-1);
       
       fai(i,j,k+1)=(-fai(i,j+1,k)*P-fai(1,j,k)*Q +5.625*R-fai(i-1,j,k+1)*S)/T; 
       
       end
       
                     
       P=(A-(1/mio)*(sai(1,j,k+1)-sai(i-1,j,k+1))/2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j-1));
       Q=(B+(1/mio)*(sai(i,j+1,k+1)+sai(i,j,k+1)*(alfa^2-1)-sai(i,j-1,k+1)*alfa^2)/2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j));
       R=((2/(alfa+1)*matrix_dr(1,j-1)^2)-alfa/(matrix_r(1,j)*(alfa+1)*matrix_dr(1,j-1))-(1/mio)*(-alfa*(sai(1,j,k+1)-sai(i-1,j,k+1)))/2*matrix_r(1,j)*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j-1));
       S=B-(1/mio)*(sai(i,j+1,k+1)+(alfa^2-1)*sai(i,j,k+1)-(alfa^2)*sai(i,j-1,k+1))/(2*matrix_r(1,j)*alfa*(alfa+1)*matrix_dtheta(1,i)*matrix_dr(1,j));
       T=D-(1/mio)*(alfa-1)*(sai(1,j,k+1)-sai(i-1,j,k+1))/2*matrix_r(1,j)*(alfa)*matrix_dtheta(1,i)*matrix_dr(1,j-1);
       
       fai(i,j,k+1)=(-fai(i,j+1,k)*P-fai(1,j,k+1)*Q-fai(i,j-1,k+1)*R-fai(i-1,j,k+1)*S)/T; 
                       
  
    
         end
  
    end
    
end




e(:,:,k+1)=abs(sai(:,:,k+1)-sai(:,:,k));
%  k=k+1;
end
%%
[x,y]=meshgrid(matrix_r,matrix_theta);
[b,o]=pol2cart(y,x);
contour(b,o,sai(:,:,k_max),100)
%%
figure
contour(b,o,fai(:,:,k_max),100)
%%
plot(matrix_r,sai(:,:,k_max));
figure
plot(matrix_r,fai(:,:,k_max));

end






