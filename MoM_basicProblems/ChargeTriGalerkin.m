clear;close;clc;

wireLength = 1;
wireRadius = 0.001;
NoOfSegments= 100;
eps0 = 8.8542e-12;
delta = wireLength/NoOfSegments;
subPoints = 3;

basisPoints = 0:delta:wireLength;

tempArray = linspace(-delta,delta,subPoints);
triangularBasis = 1-abs(tempArray)/delta;

Z_withBasisInt = zeros;
Zmn = zeros;
for m = 2:length(basisPoints)-1
    tempArray_testing = linspace(basisPoints(m-1),basisPoints(m+1),subPoints);
    for n = 2:length(basisPoints)-1
        tempArray_basis = linspace(basisPoints(n-1),basisPoints(n+1),subPoints);
        for p = 1:length(tempArray_testing)
            dist = sqrt((tempArray_basis(p)-tempArray_testing).^2 +wireRadius^2);
            Z_withBasis = triangularBasis./dist;
            Z_withBasisInt(n,p) = simp(Z_withBasis.', tempArray_basis(2)-tempArray_basis(1));
        end
        Z_withTesting = triangularBasis.*Z_withBasisInt(n,:);
        Zmn(m,n) = simp(Z_withTesting.',tempArray_testing(2)-tempArray_testing(1));
    end
end

Zmn(1,:) = []; Zmn(:,1) = [];

Vm = zeros;
Voltage = 4*pi*eps0;
for m = 1:length(basisPoints)-2
    Vm_withTesting = triangularBasis.*Voltage;
    Vm(m) = simp(Vm_withTesting.',tempArray_testing(2)-tempArray_testing(1));
end

Vm = Vm';
alpha = Zmn\Vm;
plot(basisPoints(2:length(basisPoints)-1), alpha,'r','LineWidth',2);


function q = simp(y,dx)
N = length(y);
mul = ones(1,N);
mul(2:2:N-1) = 4;
mul(3:2:N-2) = 2;

q = (dx/3)*sum(y.'.*mul);
end
