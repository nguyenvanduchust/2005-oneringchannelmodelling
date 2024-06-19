%optimaliserer phi i 2-dim.
%--------------------------------------
header;
tic;

N=80;
fD=91;
tau=0;
kappa=0;
K=0;
alphaBS=2*pi*90/360;
alphaMS=2*pi*90/360;
phimax=2*pi*2/360;
mu=2*pi*180/360;
phimax2=0.1;
a1=-pi;
p=2;

r=rand(1,N);
phi=zeros(N,1);
for n=1:N
    stop=1;
    i=0;
    sum=0;
    while stop~=0
        a2=a1+i*phimax2;
        sum=sum+phimax2*Mises(a2,kappa,mu);
       if sum>=r(n)
           phi(n,1)=a1+i*phimax2;
            stop=0;
end
        i=i+1;
end
end
%load('2dim_data.mat');

options = optimset('Display','iter','TolFun',1*10^(-12),'MaxIter',1000,'MaxFunEvals',100000000,'TypicalX',phi,'TolX',1*10^(-12));
x=fminsearch('errorfunc_2dim',phi,options);   
phi=x;
savefile='2dim_data.mat';
save(savefile,'phi');
phi_r=linspace(-pi,pi,100);
deltaBS=linspace(0,30,100);
deltaMS=linspace(0,3,100);
rr=zeros(length(deltaBS),length(deltaMS));
ss=rr;
E=rr;
for i=1:length(deltaMS)
    for j=1:length(deltaBS)
        a=exp(sqrt(-1).*pi.*deltaBS(j).*(cos(alphaBS)+phimax.*sin(alphaBS).*sin(phi_r)));
        b=exp(sqrt(-1).*pi.*deltaMS(i).*cos(phi_r-alphaMS));
        rr(i,j)=trapz(phi_r,a.^2.*b.^2.*Mises(phi_r,kappa,mu));
        sum=0;
        for n=1:N
            aa=exp(sqrt(-1).*pi.*deltaBS(j).*(cos(alphaBS)+phimax.*sin(alphaBS).*sin(phi(n))));
            bb=exp(sqrt(-1).*pi.*deltaMS(i).*cos(phi(n)-alphaMS));
            sum=sum+aa.^2.*bb.^2;
        end
        ss(i,j)=sum/N;
        E(i,j)=real(abs(rr(i,j)-ss(i,j)));
    end
end
surf(deltaBS,deltaMS,real(E))
set(gca,'YDir','reverse');






  