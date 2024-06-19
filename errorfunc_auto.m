function [E]=errorfunc_2dim(phi)
header;
taumax=0.08;
tau=linspace(0,taumax,50);
phi2=linspace(-pi,pi,100);
for i=1:length(tau)
    rr=trapz(phi2,exp(-sqrt(-1).*2.*pi.*fD.*cos(phi2-alphav).*tau(i)).*Mises(phi2,kappa,mu));
    sum=0;
    for n=1:N
        sum=sum+exp(-sqrt(-1)*2*pi*fD*cos(phi(n)-alphav)*tau(i));
    end
    ss=sum/N;
    int1(i)=(abs(rr-ss))^p;    
end
E=(1/taumax*trapz(tau,int1))^(1/p);