clc;clear;close all;
wireLength = 1;
noOfSeg = 100;
noOfBasisPoints = noOfSeg+1;
basisPoint = linspace(0,1,noOfBasisPoints);
a = 0.001;      %wire radius
delta = basisPoint(2)-basisPoint(1);
eps0 = 8.854e-12;

Z = zeros;
for m = 2:length(basisPoint)-1
    for n = 2:length(basisPoint)-1
        tempArray = linspace(basisPoint(n)-delta/2,basisPoint(n)+delta/2,5);
        Z(m,n) = simp((1./sqrt(a^2 + (basisPoint(m)-tempArray).^2)).',tempArray(2)-tempArray(1));
    end
end
Z(1,:)=[];Z(:,1)=[];

coeff = Z\(4*pi*eps0*ones(length(Z),1));

plot(basisPoint(2:length(basisPoint)-1),coeff*1e12,'-or','LineWidth',1.4);
ylim([7 15]);
xlabel('Length(m)','FontSize',12);
ylabel('Charge density (pC/m)','FontSize',12);
title('Charge density of a thin wire','FontSize',14);

function q=simp(y,dx)
N=length(y);

mul=ones(1,N);
mul(2:2:N-1)=4;
mul(3:2:N-2)=2;

q=(dx/3)*sum(y.'.*mul);
end