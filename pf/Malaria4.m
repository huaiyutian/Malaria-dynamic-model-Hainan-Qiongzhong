function ydot = Malaria4(t,y,theta,xdata)

ei1=1/theta(1); %Time from exposed to infected
i1i2=1/theta(2);%Time from symptomatic to asymptomatic
s2s1=1/theta(3);%Time of recovery from subpatient infection to susceptible
sf=theta(4);%Infectivity of subpatent infected people
ps=theta(5); %reporting rate
si2=theta(6);%Superinfection fraction, from S2 to I1
ts=theta(7);%Treatment success

s=theta(8)+theta(9)+theta(10)+theta(11)+theta(12)+theta(13);
s1=0.01*theta(8)/s;
s2=0.01*theta(9)/s;
s3=0.01*theta(10)/s;
s4=0.01*theta(11)/s;
s5=0.01*theta(12)/s;
s6=0.01*theta(13)/s;

tau=1/11;
beta_nino=theta(14);
spray=theta(15);
net=theta(16);
beta1=spline([2:2:12],[s1,s2,s3,s4,s5,s6],[1:1:12]);
beta2=repmat(beta1',16,1);


N(1)=193000;
HI1(1)=N(1)*y(1);
HE(1)=HI1(1);
HS2(1)=N(1)*y(2);
HS1(1)=N(1)-HI1(1)-HE(1)-HS2(1);
HR(1)=0;

k(1)=0;
lambda(1)=0;

birth=0.01584/360;%human population
death=0.03514/360;%human population

DrugA=repelem(xdata(ceil(t),14)/30,30);%durg-all
DrugT=repelem(xdata(ceil(t),15)/30,30);%durg-target
Insecticide=xdata(ceil(t),22);%insecticide
Net=xdata(ceil(t),25);%Mosquito net
nino3=xdata(ceil(t),20);

Cli=beta_nino*nino3;
beta3=(beta2).*(Cli-spray*Insecticide-net*Net);
beta3(beta3<0)=0;
beta=repelem(beta3,30);

j=1;
n=1;
delta=1;
for i = 2:5760 %

HS1(i)=HS1(i-1)+(N(i-1)*birth-lambda(i-1)*HS1(i-1)+HS2(i-1)*s2s1-DrugA(i-1)+(1/30)*HR(i-1)-HS1(i-1)*death); %+HI1(i-1)*ts

HE(i)=HE(i-1)+(lambda(i-1)*HS1(i-1)-ei1*HE(i-1)-HE(i-1)*death);  
if HE(i) < 0
    HE(i)=0;
end
HI1(i)=HI1(i-1)+(ei1*HE(i-1)-HI1(i-1)*ps*ts*(1/8)-HI1(i-1)*(ps)*(1-ts)*i1i2-HI1(i-1)*(1-ps)*i1i2+si2*lambda(i-1)*HS2(i-1)-HI1(i-1)*death);%si1*lambda(i-1)*HI2(i-1)+
if HI1(i) < 0
    HI1(i)=0;
end
HS2(i)=HS2(i-1)+HI1(i-1)*(ps)*(1-ts)*i1i2+HI1(i-1)*(1-ps)*i1i2-HS2(i-1)*s2s1-si2*lambda(i-1)*HS2(i-1)-DrugT(i-1)-HS2(i-1)*death;%+HI2(i-1)*i2s2 -si3*lambda(i-1)*HS2(i-1)
if HS2(i) < 0
    HS2(i)=0;
end
HR(i)=HR(i-1)+(DrugA(i-1)+DrugT(i-1)-(1/30)*HR(i-1)-HR(i-1)*death +HI1(i-1)*ps*ts*(1/8) );
if HR(i) < 0
    HR(i)=0;
end
N(i)=HS1(i-1)+HE(i-1)+HI1(i-1)+HS2(i-1)+HR(i-1);)

k(i)=k(i-1)+2*beta(i-1)*(HI1(i-1)+sf*HS2(i-1))*tau/N(i-1) - 2*tau*k(i-1);
if k(i)<0;
    k(i)=0;
end
lambda(i)=lambda(i-1)+2*tau*k(i-1)-2*tau*lambda(i-1);
if lambda(i)<0;
    lambda(i)=0;
end

reportedHI(i)=ei1*ps*HE(i-1);
allinfection(i)=lambda(i-1)*HS1(i-1);

if mod(i,30)==0   
HIObs(j)=sum(reportedHI(i-29:i));
NN(j)=N(i);
HSS1(j)=HS1(i);
HRR(j)=HR(i);
j=j+1;
end

if mod(i,360)==0
HI(n)=sum(allinfection(i-359:i-30));
n=n+1;
end
end   
HI(17:192)=0;

ydot=[HIObs(:),HI(:),beta3(:)];
