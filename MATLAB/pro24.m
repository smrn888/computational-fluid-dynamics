clear all
clc
format long
delta_x0 = input('Please Enter Delta X : ');
delta_t = input('Please Enter Delta T: ');
tf= input('Please Enter the calculation time: ');
% tf=0.5;
% delta_x0=0.01;
% delta_t=0.0001;

r=input('uniform mesh(0)/ non uniform(1) : ');
% r=1;
if r==0
    %uniform
    %delta_x = input('Please Enter Delta X : ');
    Imax=1/delta_x0;
    h=zeros(1,Imax+1);


    for i=1:Imax
        h(1)=0;
        h(:,i+1)=h(i)+delta_x0;
    end
    u=sin(pi*h);
    v=pi*cos(pi*h);
    delta_x=ones(Imax)*delta_x0;
end
if r==1
    %%   nonuniform
    %delta_x = input('Please Enter Delta X : ');
    alpha = input('Please Enter alpha : ');
    alpha=1.01;
    n=1;
    m=1;
    h1(1)=0.5;
    h2(1)=0.5;
    while h1(n)>0
        h1(n+1)=h1(n)-delta_x0/(alpha^(n-1));
        n=n+1;
    end
    while h2(m)<1
        h2(m+1)=h2(m)+delta_x0/(alpha^(m-1));
        m=m+1;
    end
    H1=zeros(n,1);
    for i=1:n-1
        H1(i)=h1(i);
    end
    H2=zeros(m-1,1);
    for i=1:m-1

        H2(i)=h2(i+1);
    end
    H2(m-1)=1;
        M1=zeros(n-1,1);
    for i=2:n-1
        M1(i)=h1(i);
    end
    M2=zeros(m-1,1);
        for i=1:m-1

        M2(i)=h2(i+1);
    end
    M2(m-1)=1;
    H=sort(cat(1,H1,H2));
    M=sort(cat(1,M1,M2));
    Imax=size(H,1)-1;
    u=sin(pi*H);
    v=pi*cos(pi*H);

    delta_x=abs(M-0.5);
         

    
end
%update
%% coeffecient

A=zeros(2,2,Imax);
B=zeros(2,2,Imax+1);
C=zeros(2,2,Imax);
R=zeros(2,1,Imax+1);
D=zeros(2,1,Imax+1);
% norm=ones(Imax+1,1);

Iteration=tf/delta_t;
% (max(norm)>(10^-6))
  for i=1:Iteration+1
% for i=1:iteration
    for i=1:Imax-1
        B(:,:,i+1)=[+1/(2*delta_t)+v(i+1)/4,-1/(2*delta_x(i))+u(i+1)/4;-1,-delta_x(i)/2];
    end
    for i=1:Imax
        A(:,:,i+1)=[(+1/(2*delta_t))+v(i)/4,+(1/(2*delta_x(i)))+(u(i)/4);0,0];
    end
    
    for i=0:Imax-1
        C(:,:,i+1)=[0,0;1,-delta_x(i+1)/2];
    end
    for i=1:Imax-1
        R(:,:,i+1)=[(-u(i)*v(i)/2)+(-u(i+1)*v(i+1)/2)+(v(i+1)/delta_x(i))-(v(i)/delta_x(i));((delta_x(i+1)/2)*(v(i+1)+v(i+2)))+u(i+1)-u(i+2)];
    end
    R(:,:,1)=[0;(delta_x(1)/2)*(v(1)+v(2))+u(1)-u(2)];
    R(:,:,Imax+1)=[-(u(Imax)*v(Imax)/2)-(u(Imax+1)*v(Imax+1)/2)+(v(Imax+1)/delta_x(Imax))-(v(Imax)/delta_x(Imax));0];
    B(:,:,1)=[1,0;-1,-delta_x(1)/2];
    B(:,:,Imax+1)=[+1/(2*delta_t)+v(Imax+1)/4,-1/(2*delta_x(Imax))+(u(Imax+1)/4);1,0];
    
    %update
    
    %% Thomas algorithm
    for i=1:Imax
        % x=(-1)*(A(:,:,i+1)/B(:,:,i));
        B(:,:,i+1)=-A(:,:,i+1)*(inv(B(:,:,i))*C(:,:,i))+B(:,:,i+1);
%         A(:,:,i+1)=-A(:,:,i+1)*(inv(B(:,:,i))*B(:,:,i))+A(:,:,i+1);
        R(:,:,i+1)=-A(:,:,i+1)*(inv(B(:,:,i))*R(:,:,i))+R(:,:,i+1);
    end
    %update
    %%delta
    D(:,:,Imax+1)=inv(B(:,:,Imax+1))*R(:,:,Imax+1);
    for i=Imax:-1:1
        D(:,:,i)=inv(B(:,:,i))*((R(:,:,i)-C(:,:,i)*D(:,:,i+1)));
    end
    %D(1,1,1)=0;
    %D(1,1,Imax+1)=0;
    %     D(:,:,1)=inv(B(:,:,1))*(R(:,:,1)-C(:,:,1)*D(:,:,2));
    %update
    %% velocity updates
    for i=1:Imax+1
        % norm(i)=sqrt(D(1,1,i)^2+D(2,1,i)^2);
        u(i)=D(1,:,i)+u(i);
        v(i)=D(2,:,i)+v(i);
    end
    E=max(max(abs(D)))
    
    %update
end
if r==0
plot(h,u,'r');
title('Velocity Profile');
xlim([0 1]);
ylim([0 1]);
xlabel('X');
ylabel('U');
end
if r==1
plot(M,u, 'r');
title('Velocity Profile');
xlim([0 1]);
xlabel('X');
ylabel('U');
end
