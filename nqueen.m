%Solve nqueen problem.

n=8;
R=zeros(m,n);
num=1;
m=100;
while R(m,1)==0
P=(1/n)*ones(n,n);
a=zeros(2*n-1,1);
for i=1:n
    a(i)=i/n;
end
for i=n+1:2*n-1
    a(i)=(2*n-i)/n;
end
b=a;
y=1;
z=1;
Larray=zeros(1,1000);

for t=1:15000
    s=randi(n,[1,4]);
    if s(1)==s(3)
        continue;
    end
    if s(2)==s(4)
        continue;
    end
    i1=s(1);
    j1=s(2);
    i2=s(3);
    j2=s(4);
    x1=P(i1,j1);
    x2=P(i1,j2);
    x3=P(i2,j1);
    x4=P(i2,j2);
    x=[x1 x2 x3 x4];
    delta1=min(1-max(x1,x4),min(x2,x3));
    delta2=-min(min(x1,x4),1-max(x2,x3));
    if i1-j1~=i2-j2 & i1+j1~=i2+j2
        k1=i1-j1+n;
        k2=i1-j2+n;
        k3=i2-j1+n;
        k4=i2-j2+n;
        k5=i1+j1-1;
        k6=i1+j2-1;
        k7=i2+j1-1;
        k8=i2+j2-1;
        u=[a(k1),a(k2),a(k3),a(k4)];
        v=[b(k5),b(k6),b(k7),b(k8)];
        if costinc1(delta1,x,u,v)>0
%            y=y+2*(P(i1,j1)+P(i2,j2)-P(i1,j2)-P(i2,j1))*delta1+4*delta1^2; 
            z=z+costinc1(delta1,x,u,v);
            P(i1,j1)=P(i1,j1)+delta1;
            P(i2,j2)=P(i2,j2)+delta1;
            P(i1,j2)=P(i1,j2)-delta1;
            P(i2,j1)=P(i2,j1)-delta1;
            a(k1)=a(k1)+delta1;
            a(k4)=a(k4)+delta1;
            a(k2)=a(k2)-delta1;
            a(k3)=a(k3)-delta1;
            b(k5)=b(k5)+delta1;
            b(k8)=b(k8)+delta1;
            b(k6)=b(k6)-delta1;
            b(k7)=b(k7)-delta1; 
            continue;
        end
        if costinc1(delta2,x,u,v)>0
 %           y=y+2*(P(i1,j1)+P(i2,j2)-P(i1,j2)-P(i2,j1))*delta2+4*delta2^2;  
 %           z=z+costinc1(delta2,x,u,v);
            P(i1,j1)=P(i1,j1)+delta2;
            P(i2,j2)=P(i2,j2)+delta2;
            P(i1,j2)=P(i1,j2)-delta2;
            P(i2,j1)=P(i2,j1)-delta2;
            a(k1)=a(k1)+delta2;
            a(k4)=a(k4)+delta2;
            a(k2)=a(k2)-delta2;
            a(k3)=a(k3)-delta2;
            b(k5)=b(k5)+delta2;
            b(k8)=b(k8)+delta2;
            b(k6)=b(k6)-delta2;
            b(k7)=b(k7)-delta2;
        end        
    else
        if i1-j1==i2-j2
            k1=i1-j1+n;
            k2=i1-j2+n;
            k3=i2-j1+n;
            k4=i1+j2-1;
            k5=i1+j1-1;
            k6=i2+j2-1;
        end
        if i1+j1==i2+j2
            k1=i2-j1+n;
            k2=i1-j1+n;
            k3=i2-j2+n;
            k4=i1+j1-1;
            k5=i1+j2-1;
            k6=i2+j1-1;
        end
        il=min(i1,i2);
        ih=max(i1,i2);
        jl=min(j1,j2);
        jh=max(j1,j2);
        u=[a(k1) a(k2) a(k3)];
        v=[b(k4) b(k5) b(k6)];
        delta=0;
        x1=P(il,jl);
        x2=P(il,jh);
        x3=P(ih,jl);
        x4=P(ih,jh);
        x=[x1 x2 x3 x4];
        delta1=min(1-max(x1,x4),min(x2,x3));
        delta2=-min(min(x1,x4),1-max(x2,x3));
        if costinc2(delta1,x,u,v)>0
            delta=delta1;
        else
            if costinc2(delta2,x,u,v)>0
                delta=delta2;
            end
        end
        if delta~=0
            P(il,jl)=P(il,jl)+delta;
            P(ih,jh)=P(ih,jh)+delta;
            P(il,jh)=P(il,jh)-delta;
            P(ih,jl)=P(ih,jl)-delta;
            a(k1)=a(k1)+2*delta;
            a(k2)=a(k2)-delta;
            a(k3)=a(k3)-delta;
            b(k4)=b(k4)-2*delta;
            b(k5)=b(k5)+delta;
            b(k6)=b(k6)+delta;
        end
    end
%    Larray(t)=Lvalue(P,a,b,n);
%    disp(trace((P'*P)));
end
if max(a)<=1 &max(b)<=1
    result=zeros(1,n);
    [MAX,result]=max(P);
    if MAX>0.5
        R(num,:)=result;
        num=num+1;
    end
end
end

%Result statistics.
if n<10
    S=zeros(m,1);
    for i=1:m
        for j=1:n
            S(i)=S(i)+R(i,j)*10^(n-j);
        end
    end
end

function C=costinc1(delta,x,u,v)
C=(2*(x(1)+x(4)-x(2)-x(3))*delta+4*delta^2)*1;
d=1;
if d==1
C=C+max(u(1)-1,0)-max(u(1)+delta-1,0);
C=C+max(u(4)-1,0)-max(u(4)+delta-1,0);
C=C+max(u(2)-1,0)-max(u(2)-delta-1,0);
C=C+max(u(3)-1,0)-max(u(3)-delta-1,0);
C=C+max(v(1)-1,0)-max(v(1)+delta-1,0);
C=C+max(v(4)-1,0)-max(v(4)+delta-1,0);
C=C+max(v(2)-1,0)-max(v(2)-delta-1,0);
C=C+max(v(3)-1,0)-max(v(3)-delta-1,0);
end
end

function C=costinc2(delta,x,u,v)
C=(2*(x(1)+x(4)-x(2)-x(3))*delta+4*delta^2)*1;
C=C+max(u(1)-1,0)-max(u(1)+2*delta-1,0);
C=C+max(u(2)-1,0)-max(u(2)-delta-1,0);
C=C+max(u(3)-1,0)-max(u(3)-delta-1,0);
C=C+max(v(1)-1,0)-max(v(1)-2*delta-1,0);
C=C+max(v(2)-1,0)-max(v(2)+delta-1,0);
C=C+max(v(3)-1,0)-max(v(3)+delta-1,0);
end

function L=Lvalue(P,a,b,n)
L=trace(P'*P);
for i=1:2*n-1
    L=L-max(a(i)-1,0)-max(b(i)-1,0);
end
end