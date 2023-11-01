clc;clear;close all;
noOfSeg = 35;
noOfBasisPoints = noOfSeg+1;
plateLength = 1;
basisPoint = linspace(0,plateLength,noOfBasisPoints);
delta = basisPoint(2)-basisPoint(1);
basisPoint_midpoint = basisPoint(1:length(basisPoint)-1)+delta/2;
eps0 = 8.854e-12;

label = 1;
for m = 1:length(basisPoint_midpoint)
    for n = 1:length(basisPoint_midpoint)
        xy_cord(label,:) = [basisPoint_midpoint(m) basisPoint_midpoint(n)];
        label = label+1;
    end
end

for m = 1:length(xy_cord)
    for n = 1:length(xy_cord)
        if m==n
            Z(m,n) = (delta/(pi*eps0))*0.8841;
        else
            Z(m,n) = (((delta^2)/(4*pi*eps0*sqrt((xy_cord(m,1)-xy_cord(n,1))^2 +(xy_cord(m,2)-xy_cord(n,2))^2))));
        end
    end
end

coeff = Z\(ones(length(Z),1));

for m = 1:sqrt(length(coeff))
    coeff_diag(m) = coeff(noOfSeg*(m-1)+m,1);
end

plot(0:length(coeff_diag)-1,(coeff_diag)*1e12,'LineWidth',1.4);
