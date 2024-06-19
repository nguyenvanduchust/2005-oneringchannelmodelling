%optimaliserer phi i 1-dim.
%--------------------------------------
header;
tic;

N=25;
fD=91;
tau=0;
kappa=100;
K=0;
alphaBS=2*pi*90/360;
alphaMS=2*pi*90/360;
phimax=2*pi*2/360;
alphav=2*pi*180/360;
mu=2*pi*180/360;
phimax2=0.1;
a1=-pi;
p=2;

%r=rand(1,N);
%phi=linspace(1,N);
%for n=1:N
%    stop=1;
%    i=0;
%    sum=0;
%    while stop~=0
%        a2=a1+i*phimax2;
%        sum=sum+phimax2*Mises(a2,kappa,mu);
%        if sum>=r(n)
%            phi(n,1)=a1+i*phimax2;
%            stop=0;
%        end
%        i=i+1;
%    end
%end
load('auto_data.mat');

options = optimset('Display','iter','TolFun',1*10^(-20),'MaxIter',10000,'MaxFunEvals',100000000,'TypicalX',phi,'TolX',1*10^(-20));
x=fminsearch('errorfunc_auto',phi,options);   
phi=x;
savefile='auto_data.mat';
save(savefile,'phi');
tau=linspace(0,0.1,100);
phi2=linspace(-pi,pi,200);
for i=1:length(tau)
    rr(i)=trapz(phi2,exp(-sqrt(-1)*2*pi*fD.*cos(phi2-alphav)*tau(i)).*Mises(phi2,kappa,mu));
    sum=0;
    for n=1:N
        sum=sum+exp(-sqrt(-1)*2*pi*fD*cos(phi(n)-alphav)*tau(i));
    end
    ss(i)=sum/N;   
end
plot(tau,real(rr),tau,real(ss))





  