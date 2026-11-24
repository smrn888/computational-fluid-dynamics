clear all
clc
format long
d_r=input('delta r : ');
d_teta=input('delta teta : ');
alfa=input('alfa : ');
i_max=floor(2*pi/d_teta);
j=0;
y=0.3;
while(0.5>y)
    y=y+d_r*alfa^(j+1);
    j=j+1;
end
j_max=j+1;
r=zeros(j_max,1);
dr=zeros(j_max-1,1);
dr(1)=d_r;
teta=zeros(i_max,1);
r(1,1)=.3;
for i=1:j_max-1
    r(i+1,1)=r(i,1)+d_r*alfa^(i-1);
end
for i=1:j_max-1
    dr(i)=d_r*alfa^(i-1);
end
r(j_max)=.5;
dr(j_max-1)=.5-r(j_max-1);
%%
for i=1:i_max
    teta(i)=i*d_teta;
end
%%
say_0=zeros(i_max,1);
say=zeros(i_max,j_max);
vort=ones(i_max,j_max);
say_p=zeros(i_max,j_max);
vort_p=zeros(i_max,j_max);
A=zeros(j_max,1);
B=zeros(j_max,1);
C=zeros(j_max,1);
D=zeros(j_max,1);
H=zeros(i_max,j_max);
w=5;
   v=.1;
   e=1;
   k=0;
 vort(:,j_max)=(-5.625);
 vort(:,1)=(-5.625);
say(:,j_max)=1;
X=zeros(i_max,j_max);
%%
% for k=1:1000
 while(e>10^-6)
     k=k+1;
       say_p=say;
       vort_p=vort;
 for j=2:j_max
        A(j,1)=-1/[r(j)*((alfa^2)+alfa)*dr(j-1)]-2/[((alfa^2)+alfa)*dr(j-1)];
        B(j,1)=(alfa^2)/[r(j)*(alfa^2+alfa)*dr(j-1)]-(2*alfa)/[(alfa^2+alfa)*(dr(j-1))^2];
        C(j,1)=-1/[r(j)^2*d_teta^2];
        D(j,1)=-[1-alfa^2]/[r(j)*(alfa^2+alfa)*dr(j-1)]-2*(1+alfa)/[(alfa^2+alfa)*(dr(j-1)^2)]-2/[r(j)^2*d_teta^2];
 end
for i=1:i_max
     if i==1
         say(i,2)=(-vort(i,2)+say(i,3)*A(2)+(d_r*0.3*5)*B(2)+say(i_max,2)*C(2)+say(i+1,2)*C(2))/(D(2)+B(2));
     elseif i~=i_max
         say(i,2)=(-vort(i,2)+say(i,3)*A(2)+(d_r*0.3*5)*B(2)+say(i-1,2)*C(2)+say(i+1,2)*C(2))/(D(2)+B(2));
     else
         say(i,2)=(-vort(i,2)+say(i,3)*A(2)+(d_r*0.3*5)*B(2)+say(i-1,2)*C(2)+say_p(1,2)*C(2))/(D(2)+B(2));
     end
end
for j=3:j_max-1
    for i=1:i_max
        if i==1
            say(1,j)=(-vort(1,j)+say(1,j+1)*A(j)+say(1,j-1)*B(j)+say(i_max,j)*C(j)+say(2,j)*C(j))/(D(j));
        elseif i~=i_max
            say(i,j)=(-vort(i,j)+say(i,j+1)*A(j)+say(i,j-1)*B(j)+say(i-1,j)*C(j)+say(i+1,j)*C(j))/(D(j));
        else
            say(i_max,j)=(-vort(i_max,j)+say(i_max,j+1)*A(j)+say(i_max,j-1)*B(j)+say(i_max-1,j)*C(j)+say_p(1,j)*C(j))/(D(j));
        end
    end
end          
for i=1:i_max
   say(i,1)=say(i,2)+r(1)*5*dr(1);
end      
for i=1:i_max
   say(i,j_max)=1;
end   

for i=1:i_max
    for j=2:j_max-1
        if i==1
 
    X(i,j)=(say(i+1,j)-say(i_max,j))/(2*r(j)*d_teta*dr(j-1)*(alfa+1))-v/(r(j)*dr(j-1)*(alfa+1))-...
        v/[(dr(j-1))^2*alfa*(alfa+1)];
    Y(i,j)=[say(i_max,j)-say(i+1,j)]/[2*r(j)*d_teta*dr(j-1)*(alfa+1)]+v/[r(j)*dr(j-1)*(alfa+1)]-...
        v/[(dr(j-1))^2*(alfa+1)];
    Z(i,j)=[say(i,j-1)-say(i,j+1)]/[2*r(j)*d_teta*dr(j-1)*(alfa+1)]-v/([r(j)*d_teta]^2);
    T(i,j)=[say(i,j+1)-say(i,j-1)]/[2*r(j)*d_teta*dr(j-1)*(alfa+1)]-v/([r(j)*d_teta]^2);
    S(i,j)=-v/[(dr(j-1))^2*alfa]-2*v/([r(j)*d_teta]^2);
    vort(i,j)=[X(i,j)*vort(i,j+1)+Y(i,j)*vort(i,j-1)+Z(i,j)*vort(i+1,j)+T(i,j)*vort(i_max,j)]/S(i,j);
    %%
        elseif i~=i_max
    X(i,j)=(say(i+1,j)-say(i-1,j))/(2*r(j)*d_teta*dr(j-1)*(alfa+1))-v/(r(j)*dr(j-1)*(alfa+1))-...
        v/[(dr(j-1))^2*alfa*(alfa+1)];
    Y(i,j)=[say(i-1,j)-say(i+1,j)]/[2*r(j)*d_teta*dr(j-1)*(alfa+1)]+v/[r(j)*dr(j-1)*(alfa+1)]-...
        v/[(dr(j-1))^2*(alfa+1)];
    Z(i,j)=[say(i,j-1)-say(i,j+1)]/[2*r(j)*d_teta*dr(j-1)*(alfa+1)]-v/([r(j)*d_teta]^2);
    T(i,j)=[say(i,j+1)-say(i,j-1)]/[2*r(j)*d_teta*dr(j-1)*(alfa+1)]-v/([r(j)*d_teta]^2);
    S(i,j)=-v/[(dr(j-1))^2*alfa]-2*v/([r(j)*d_teta]^2);
    vort(i,j)=[X(i,j)*vort(i,j+1)+Y(i,j)*vort(i,j-1)+Z(i,j)*vort(i+1,j)+T(i,j)*vort(i-1,j)]/S(i,j);
    %%
        else

    X(i,j)=(say(1,j)-say(i-1,j))/(2*r(j)*d_teta*dr(j-1)*(alfa+1))-v/(r(j)*dr(j-1)*(alfa+1))-...
        v/[(dr(j-1))^2*alfa*(alfa+1)];
    Y(i,j)=[say(i-1,j)-say(1,j)]/[2*r(j)*d_teta*dr(j-1)*(alfa+1)]+v/[r(j)*dr(j-1)*(alfa+1)]-...
        v/[(dr(j-1))^2*(alfa+1)];
    Z(i,j)=[say(i,j-1)-say(i,j+1)]/[2*r(j)*d_teta*dr(j-1)*(alfa+1)]-v/([r(j)*d_teta]^2);
    T(i,j)=[say(i,j+1)-say(i,j-1)]/[2*r(j)*d_teta*dr(j-1)*(alfa+1)]-v/([r(j)*d_teta]^2);
    S(i,j)=-v/[(dr(j-1))^2*alfa]-2*v/([r(j)*d_teta]^2);
    vort(i,j)=[X(i,j)*vort(i,j+1)+Y(i,j)*vort(i,j-1)+Z(i,j)*vort(1,j)+T(i,j)*vort(i-1,j)]/S(i,j);
        end
        end
end
say
vort
 say-say_p;
vort-vort_p;
e_say=say-say_p;
e_vort=vort-vort_p;
e1=max(max(abs(e_say)));
e2=max(max(abs(e_vort)));
e(1,k)=max(e1,e2);
 end
  [x,y]=meshgrid(r,teta);
[b,o]=pol2cart(y,x);
contour(b,o,say,100)
%%
figure
contour(b,o,vort,100)
%%
plot(r,say);
figure
plot(r,vort);

e;
   