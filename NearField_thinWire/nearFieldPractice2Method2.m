clear;close;clc;tic;

noOfSegments = 32;
noOfBasisPoints = noOfSegments+1;
lambda = 1;
wireLength = lambda/2;
wireRadius = 0.005*lambda;
k = 2*pi/lambda;
eps0 = 8.8542e-12;
mu0 = pi*4e-7;
c0 = 1/sqrt(mu0*eps0);
f0 = c0/lambda;
omega = 2*pi*f0;
delta = wireLength/noOfSegments;
subPoints = 5;

Vm = zeros((noOfSegments-1),1);
Vm((noOfSegments/2),1) = 1;

basisPoints = -wireLength/2:delta:wireLength/2;

s = linspace(0,delta,subPoints); delta_s = s(2)-s(1);
basisPoints_nMinus = zeros;
basisPoints_mPlus = zeros;
basisPoints_mMinus = zeros;
Zmn = zeros;
factor1 = 1j*omega*mu0*delta/(4*pi);
factor2 = 1/(1j*omega*4*pi*eps0);

for m = 2:noOfBasisPoints-1
    for n = 2:noOfBasisPoints-1

        basisPoints_nMinus(n) = (basisPoints(n) + basisPoints(n-1))/2;
        basisPoints_mPlus(m) = (basisPoints(m) + basisPoints(m+1))/2;
        basisPoints_mMinus(m) = (basisPoints(m) + basisPoints(m-1))/2;
        psi_n_m = calculatePsi(0,0,basisPoints,basisPoints_nMinus,s,wireRadius,k,delta,delta_s,m,n,NaN);
        psi_nPlus_mPlus = calculatePsi(0,0,basisPoints_mPlus,basisPoints,s,wireRadius,k,delta,delta_s,m,n,NaN);
        psi_nPlus_mMinus = calculatePsi(0,0,basisPoints_mMinus,basisPoints,s,wireRadius,k,delta,delta_s,m,n,NaN);
        psi_nMinus_mPlus = calculatePsi(0,0,basisPoints_mPlus,basisPoints,s,wireRadius,k,delta,delta_s,m,n-1,NaN);
        psi_nMinus_mMinus = calculatePsi(0,0,basisPoints_mMinus,basisPoints,s,wireRadius,k,delta,delta_s,m,n-1,NaN);

        Zmn(m,n) = factor1*delta*psi_n_m + factor2*(psi_nPlus_mPlus - psi_nMinus_mPlus - psi_nPlus_mMinus + psi_nMinus_mMinus);
    end
end
Zmn(1,:) = []; Zmn(:,1) = [];

I = Zmn\Vm;
I1 = I.'*1e3;
magI = abs(I1);

fig1 = figure('Color','w');
% set(fig1,'FinalDraft','off');
yyaxis left
plot(basisPoints,[0 abs(I1) 0],'LineWidth',2);
ylim([0 12]);
ylabel('Magnitude of Current (mA)','FontSize',14);
hold on;
yyaxis right
plot(basisPoints,[rad2deg(angle(I1(1))) rad2deg(angle(I1)) rad2deg(angle(I1(noOfSegments-1)))],'--','LineWidth',2);
ylim([-180 180])
ylabel('Phase of Current (Degree)','FontSize',14)
hold on;
xlim([-wireLength/2 wireLength/2])
xlabel('Wire Length (m)','FontSize',14)
legend('Magnitude','Phase','Location','northwest','FontSize',14);


z_vector = 0:0.001*lambda:0.3*lambda;
testingDipoleHalfLength = 0.001*lambda/2;
rho_vector = [0.02 0.03 0.05 0.1]*lambda;
efield_Z = zeros; angle_efield_Z = zeros;

for q = 1:length(rho_vector)
    for p = 1:length(z_vector)

        basisPoints(noOfBasisPoints+1) = z_vector(p) - testingDipoleHalfLength;
        basisPoints(noOfBasisPoints+2) = z_vector(p);
        basisPoints(noOfBasisPoints+3) = z_vector(p) + testingDipoleHalfLength;
        testingDipoleUpperMidPoint = (basisPoints(noOfBasisPoints+2)+basisPoints(noOfBasisPoints+3))/2;
        testingDipoleLowerMidPoint = (basisPoints(noOfBasisPoints+2)+basisPoints(noOfBasisPoints+1))/2;
        diffz = testingDipoleUpperMidPoint - testingDipoleLowerMidPoint;
        esum = 0;
        for n = 2:noOfBasisPoints-1
            psi_n_m = calculatePsi(0,rho_vector,z_vector,basisPoints_nMinus,s,wireRadius,k,delta,delta_s,p,n,q);
            psi_nPlus_mPlus = calculatePsi(0,rho_vector,testingDipoleUpperMidPoint,basisPoints,s,wireRadius,k,delta,delta_s,NaN,n,q);
            psi_nPlus_mMinus = calculatePsi(0,rho_vector,testingDipoleLowerMidPoint,basisPoints,s,wireRadius,k,delta,delta_s,NaN,n,q);
            psi_nMinus_mPlus = calculatePsi(0,rho_vector,testingDipoleUpperMidPoint,basisPoints,s,wireRadius,k,delta,delta_s,NaN,n-1,q);
            psi_nMinus_mMinus = calculatePsi(0,rho_vector,testingDipoleLowerMidPoint,basisPoints,s,wireRadius,k,delta,delta_s,NaN,n-1,q);

            matrix = factor1*diffz*psi_n_m + factor2*(psi_nPlus_mPlus-psi_nMinus_mPlus-psi_nPlus_mMinus+psi_nMinus_mMinus);
            esum = esum + matrix*I(n-1);
        end
        efield_Z(q,p) = abs(-esum/0.001/1.414);
        angle_efield_Z(q,p) = rad2deg(angle(-esum));
    end
end
fig2 = figure('Color','w');

plot(z_vector,efield_Z(1,:),'r','LineWidth',2);hold on;
plot(z_vector,efield_Z(2,:),'g','LineWidth',2);
plot(z_vector,efield_Z(3,:),'b','LineWidth',2);
plot(z_vector,efield_Z(4,:),'m','LineWidth',2);hold off;

%title('z-comp. of Electric Field for 0.5\lambda dipole','Fontsize',10);
xlabel('z/\lambda','FontSize',14);
ylabel('E_z in V/m','Fontsize',14);
legend('\rho/\lambda=0.02','\rho/\lambda=0.03','\rho/\lambda=0.05','\rho/\lambda=0.1','Location','north','Orientation','vertical','Fontsize',14);

fig3 = figure('Color','w');
plot(z_vector,angle_efield_Z(1,:),'-r','LineWidth',2);hold on;
plot(z_vector,angle_efield_Z(2,:),'-g','LineWidth',2);
plot(z_vector,angle_efield_Z(3,:),'-b','LineWidth',2);
plot(z_vector,angle_efield_Z(4,:),'-m','LineWidth',2);hold off;
xlabel('z/\lambda','FontSize',14);
ylabel('Phase of E_Z (in degrees)','Fontsize',14);
ylim([-180 180]);
legend('\rho/\lambda=0.02','\rho/\lambda=0.03','\rho/\lambda=0.05','\rho/\lambda=0.1','Location','north','Orientation','vertical','Fontsize',14);

rho_vector = [0.01 0.02 0.05]*lambda;
efield_Y = zeros; angle_efield_Y = zeros; efield_X = zeros;

for q = 1:length(rho_vector)
    for p = 1:length(z_vector)
        basisPoints(noOfBasisPoints+1) = z_vector(p);
        basisPoints(noOfBasisPoints+2) = z_vector(p);
        basisPoints(noOfBasisPoints+3) = z_vector(p);
        testingDipoleLowerPoint_y = rho_vector(q) - testingDipoleHalfLength;
        testingDipoleUpperPoint_y = rho_vector(q) + testingDipoleHalfLength;
        testingDipoleMainPoint_y = rho_vector(q);
        testingDipoleLowerMidPoint = (testingDipoleLowerPoint_y + testingDipoleMainPoint_y)/2;
        testingDipoleUpperMidPoint = (testingDipoleUpperPoint_y + testingDipoleMainPoint_y)/2;

        esum = 0;
        for n = 2:noOfBasisPoints-1
            psi_nPlus_mPlus = calculatePsi(0,testingDipoleUpperMidPoint,z_vector,basisPoints,s,wireRadius,k,delta,delta_s,p,n,NaN);
            psi_nPlus_mMinus = calculatePsi(0,testingDipoleLowerMidPoint,z_vector,basisPoints,s,wireRadius,k,delta,delta_s,p,n,NaN);
            psi_nMinus_mPlus = calculatePsi(0,testingDipoleUpperMidPoint,z_vector,basisPoints,s,wireRadius,k,delta,delta_s,p,n-1,NaN);
            psi_nMinus_mMinus = calculatePsi(0,testingDipoleLowerMidPoint,z_vector,basisPoints,s,wireRadius,k,delta,delta_s,p,n-1,NaN);

            matrix = factor2*(psi_nPlus_mPlus-psi_nMinus_mPlus-psi_nPlus_mMinus+psi_nMinus_mMinus);
            esum = esum + matrix*I(n-1);
        end
        efield_Y(q,p) = abs(-esum/0.001/1.414);
        angle_efield_Y(q,p) = rad2deg(angle(-esum));
    end
end

for q = 1:length(rho_vector)
    for p = 1:length(z_vector)
        basisPoints(noOfBasisPoints+1) = z_vector(p);
        basisPoints(noOfBasisPoints+2) = z_vector(p);
        basisPoints(noOfBasisPoints+3) = z_vector(p);
        testingDipoleLowerPoint_x = -testingDipoleHalfLength;
        testingDipoleUpperPoint_x = testingDipoleHalfLength;
        testingDipoleMainPoint_x = 0;
        testingDipoleLowerMidPoint = (testingDipoleMainPoint_x + testingDipoleLowerPoint_x)/2;
        testingDipoleUpperMidPoint = (testingDipoleUpperPoint_x + testingDipoleMainPoint_x)/2;

        esum = 0;
        for n = 2:noOfBasisPoints-1
            psi_nPlus_mPlus = calculatePsi(testingDipoleUpperMidPoint,rho_vector,z_vector,basisPoints,s,wireRadius,k,delta,delta_s,p,n,q);
            psi_nPlus_mMinus = calculatePsi(testingDipoleLowerMidPoint,rho_vector,z_vector,basisPoints,s,wireRadius,k,delta,delta_s,p,n,q);
            psi_nMinus_mPlus = calculatePsi(testingDipoleUpperMidPoint,rho_vector,z_vector,basisPoints,s,wireRadius,k,delta,delta_s,p,n-1,q);
            psi_nMinus_mMinus = calculatePsi(testingDipoleLowerMidPoint,rho_vector,z_vector,basisPoints,s,wireRadius,k,delta,delta_s,p,n-1,q);

            matrix = factor2*(psi_nPlus_mPlus-psi_nMinus_mPlus-psi_nPlus_mMinus+psi_nMinus_mMinus);
            esum = esum + matrix*I(n-1);
        end
        efield_X(q,p) = abs(-esum/0.001/1.414);
    end
end
efield_rho = zeros(q,p);
for q = 1:length(rho_vector)
    efield_rho(q,:) = sqrt(efield_X(q,:).^2 + efield_Y(q,:).^2);
end

fig4 = figure('Color','w');
plot(z_vector,efield_rho(1,:),'r','LineWidth',2);hold on;
plot(z_vector,efield_rho(2,:),'b','LineWidth',2);
plot(z_vector,efield_rho(3,:),'m','LineWidth',2);hold off;


xlabel('z/\lambda','Fontsize',14);
ylabel('E_\rho in V/m','Fontsize',14);
legend('\rho/\lambda=0.01','\rho/\lambda=0.02','\rho/\lambda=0.05','Location','north','Orientation','vertical','Fontsize',14);

z_vector = z_vector(2:length(z_vector));
fig5 = figure('Color','w');
plot(z_vector,angle_efield_Y(1,2:length(angle_efield_Y)),'r','LineWidth',2);hold on;
plot(z_vector,angle_efield_Y(2,2:length(angle_efield_Y)),'b','LineWidth',2);
plot(z_vector,angle_efield_Y(3,2:length(angle_efield_Y)),'m','LineWidth',2);hold off;


xlabel('z/\lambda','Fontsize',14);
ylabel('Phase of E_\rho (in degrees)','Fontsize',14);
ylim([-180 180]);
legend('\rho/\lambda=0.01','\rho/\lambda=0.02','\rho/\lambda=0.05','Location','north','Orientation','vertical','Fontsize',14);

toc;

function psi = calculatePsi(x,y,z,z_prime,s,a,k,delta,delta_s,m,n,q)
if  isnan(q)
    dist = sqrt((z(m)-z_prime(n)-s).^2 + x^2 +y^2 + a^2);
elseif isnan(m)
    dist = sqrt((z-z_prime(n)-s).^2 + x^2 +y(q)^2 + a^2);
else
    dist  = sqrt((z(m)-z_prime(n)-s).^2 + x^2 +y(q)^2 + a^2);
end
F = exp(-1j*k*dist)./dist;
psi = (1/delta)*simp(F.',delta_s);
end

function q=simp(y,dx)
N=length(y);

mul=ones(1,N);
mul(2:2:N-1)=4;
mul(3:2:N-2)=2;

q=(dx/3)*sum(y.'.*mul);
end