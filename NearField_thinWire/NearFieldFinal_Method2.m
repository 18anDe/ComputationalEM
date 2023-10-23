clear;close;clc;tic;

noOfSegments = 32;
noOfBasisPoint = noOfSegments+1;
lambda = 1;
wireLength = lambda/2;
wireRadius = 0.005*lambda;
eps0 = 8.8542e-12;
mu0 = pi*4e-7;
c0 = 1/sqrt(mu0*eps0);
freq = c0/lambda;
omega = 2*pi*freq;
k = 2*pi/lambda;
delta = wireLength/noOfSegments;

basisPoint_z = -wireLength/2:delta:wireLength/2;
basisPoint_y = zeros(length(basisPoint_z));
basisPoint_x = zeros(length(basisPoint_z));
 
Vm = zeros((noOfSegments-1),1);
Vm((noOfSegments/2),1) = 1;

s = linspace(0,delta,5);delta_s = s(2)-s(1); Zmn = zeros;
factor1 = 1j*omega*mu0*delta/(4*pi);
factor2 = 1/(1j*omega*eps0*4*pi);
basisPoint_nMinus = zeros;basisPoint_mPlus = zeros;basisPoint_mMinus = zeros;
for m = 2:length(basisPoint_z)-1
    for n = 2:length(basisPoint_z)-1

        basisPoint_nMinus(n) = (basisPoint_z(n)+basisPoint_z(n-1))/2;
        F1 = exp(-1j*k*sqrt((basisPoint_z(m)-basisPoint_nMinus(n)-s).^2 + wireRadius^2))./sqrt((basisPoint_z(m)-basisPoint_nMinus(n)-s).^2 + wireRadius^2);
        psi_n_m = (1/delta)*simp(F1.',delta_s);

        basisPoint_mPlus(m) = (basisPoint_z(m)+basisPoint_z(m+1))/2;
        F2 = exp(-1j*k*sqrt((basisPoint_mPlus(m)-basisPoint_z(n)-s).^2 + wireRadius^2))./sqrt((basisPoint_mPlus(m)-basisPoint_z(n)-s).^2 + wireRadius^2);
        psi_nPlus_mPlus = (1/delta)*simp(F2.',delta_s);

        F3 = exp(-1j*k*sqrt((basisPoint_mPlus(m)-basisPoint_z(n-1)-s).^2 + wireRadius^2))./sqrt((basisPoint_mPlus(m)-basisPoint_z(n-1)-s).^2 + wireRadius^2);
        psi_nMinus_mPlus = (1/delta)*simp(F3.',delta_s);

        basisPoint_mMinus(m) = (basisPoint_z(m) + basisPoint_z(m-1))/2;
        F4 = exp(-1j*k*sqrt((basisPoint_mMinus(m)-basisPoint_z(n)-s).^2 + wireRadius^2))./sqrt((basisPoint_mMinus(m)-basisPoint_z(n)-s).^2 + wireRadius^2);
        psi_nPlus_mMinus = (1/delta)*simp(F4.',delta_s);

        F5 = exp(-1j*k*sqrt((basisPoint_mMinus(m)-basisPoint_z(n-1)-s).^2 + wireRadius^2))./sqrt((basisPoint_mMinus(m)-basisPoint_z(n-1)-s).^2 + wireRadius^2);
        psi_nMinus_mMinus = (1/delta)*simp(F5.',delta_s);

        Zmn(m,n) = factor1*delta*psi_n_m + factor2*(psi_nPlus_mPlus-psi_nMinus_mPlus-psi_nPlus_mMinus+psi_nMinus_mMinus);
    end
end
Zmn(1,:) = []; Zmn(:,1) = [];

I = Zmn\Vm;
I1 = I.'*1e3;
magI = abs(I1);

fig1 = figure('Color','w');
% set(fig1,'FinalDraft','off');
yyaxis left
plot(basisPoint_z,[0 abs(I1) 0],'LineWidth',2);
ylim([0 12])
ylabel('Magnitude of Current (mA)','FontSize',14)
hold on;
yyaxis right
plot(basisPoint_z,[rad2deg(angle(I1(1))) rad2deg(angle(I1)) rad2deg(angle(I1(noOfSegments-1)))],'--','LineWidth',2);
ylim([-180 180])
ylabel('Phase of Current (Degree)','FontSize',14)
hold on;
xlim([-wireLength/2 wireLength/2])
xlabel('Wire Length (m)','FontSize',14)
legend('Magnitude','Phase','Location','northwest','FontSize',14);


z_vector2 = 0:0.001*lambda:0.3*lambda;
testingDipoleHalfLength = 0.001*lambda/2;
rho_vector = [0.02 0.03 0.05 0.1]*lambda;
efield_Z2 = zeros;efield_Y = zeros;efield_X = zeros;

for q = 1:length(rho_vector)
    for p = 1:length(z_vector2)

        basisPoint_z(noOfBasisPoint+1) = z_vector2(p) - testingDipoleHalfLength;
        basisPoint_z(noOfBasisPoint+2) = z_vector2(p);
        basisPoint_z(noOfBasisPoint+3) = z_vector2(p) + testingDipoleHalfLength;
        testingDipoleUpperMidPoint = (basisPoint_z(noOfBasisPoint+2) + basisPoint_z(noOfBasisPoint+3))/2;
        testingDipoleLowerMidPoint = (basisPoint_z(noOfBasisPoint+2) + basisPoint_z(noOfBasisPoint+1))/2;
        diffz = testingDipoleUpperMidPoint - testingDipoleLowerMidPoint;

        esum = 0;
        for n = 2:noOfBasisPoint-1
            basisPoint_nMinus(n) = (basisPoint_z(n)+basisPoint_z(n-1))/2;
            F1 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_nMinus(n)-s).^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_nMinus(n)-s).^2 + rho_vector(q)^2 + wireRadius^2);
            psi_n_m = (1/delta)*simp(F1.',delta_s);

            F2 = exp(-1j*k*sqrt((testingDipoleUpperMidPoint-basisPoint_z(n)-s).^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((testingDipoleUpperMidPoint-basisPoint_z(n)-s).^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nPlus_mPlus = (1/delta)*simp(F2.',delta_s);

            F3 = exp(-1j*k*sqrt((testingDipoleUpperMidPoint-basisPoint_z(n-1)-s).^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((testingDipoleUpperMidPoint-basisPoint_z(n-1)-s).^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nMinus_mPlus = (1/delta)*simp(F3.',delta_s);

            F4 = exp(-1j*k*sqrt((testingDipoleLowerMidPoint-basisPoint_z(n)-s).^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((testingDipoleLowerMidPoint-basisPoint_z(n)-s).^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nPlus_mMinus = (1/delta)*simp(F4.',delta_s);

            F5 = exp(-1j*k*sqrt((testingDipoleLowerMidPoint-basisPoint_z(n-1)-s).^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((testingDipoleLowerMidPoint-basisPoint_z(n-1)-s).^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nMinus_mMinus = (1/delta)*simp(F5.',delta_s);

            matrix = factor1*diffz*psi_n_m + factor2*(psi_nPlus_mPlus-psi_nMinus_mPlus-psi_nPlus_mMinus+psi_nMinus_mMinus);
            esum = esum + matrix*I(n-1);
        end
        efield_Z2(q,p) = abs(-esum/0.001/1.414);
        angle_efield_Z2(q,p) = rad2deg(angle(-esum));
    end
end



fig2 = figure('Color','w');

plot(z_vector2,efield_Z2(1,:),'r','LineWidth',2);hold on;
plot(z_vector2,efield_Z2(2,:),'g','LineWidth',2);hold on;
plot(z_vector2,efield_Z2(3,:),'b','LineWidth',2);hold on;
plot(z_vector2,efield_Z2(4,:),'m','LineWidth',2);hold on;

%title('z-comp. of Electric Field for 0.5\lambda dipole','Fontsize',10);
xlabel('z/\lambda','FontSize',14);
ylabel('E_z in V/m','Fontsize',14);
legend('\rho/\lambda=0.02','\rho/\lambda=0.03','\rho/\lambda=0.05','\rho/\lambda=0.1','Location','north','Orientation','vertical','Fontsize',14);
%ah1=axes('position',get(gca,'position'),'visible','off');
%title(lgd1,'Method  2','Fontsize',14)

fig3 = figure('Color','w');
plot(z_vector2,angle_efield_Z2(1,:),'-r','LineWidth',2);hold on;
plot(z_vector2,angle_efield_Z2(2,:),'-g','LineWidth',2);hold on;
plot(z_vector2,angle_efield_Z2(3,:),'-b','LineWidth',2);hold on;
plot(z_vector2,angle_efield_Z2(4,:),'-m','LineWidth',2);hold on;
xlabel('z/\lambda','FontSize',14);
ylabel('Phase of E_Z (in degrees)','Fontsize',14);
ylim([-180 180]);
legend('\rho/\lambda=0.02','\rho/\lambda=0.03','\rho/\lambda=0.05','\rho/\lambda=0.1','Location','north','Orientation','vertical','Fontsize',14);

rho_vector = [0.01 0.02 0.05]*lambda;

for q = 1:length(rho_vector)
    for p = 1:length(z_vector2)
        basisPoint_z(noOfBasisPoint+1) = z_vector2(p);
        basisPoint_z(noOfBasisPoint+2) = z_vector2(p);
        basisPoint_z(noOfBasisPoint+3) = z_vector2(p);
        testingDipoleLowerPoint_y = rho_vector(q)-testingDipoleHalfLength;
        tesingDipoleUpperPoint_y = rho_vector(q)+testingDipoleHalfLength;
        testingDipoleMainPoint_y = rho_vector(q);
        testingDipoleLowerMidPoint = (testingDipoleLowerPoint_y+testingDipoleMainPoint_y)/2;
        testingDipoleUpperMidPoint = (testingDipoleMainPoint_y+tesingDipoleUpperPoint_y)/2;
        
        esum = 0;
        for n = 2:noOfBasisPoint-1
            F2 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleUpperMidPoint^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleUpperMidPoint^2 + wireRadius^2);
            psi_nPlus_mPlus = (1/delta)*simp(F2.',delta_s);

            F3 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleUpperMidPoint^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleUpperMidPoint^2 + wireRadius^2);
            psi_nMinus_mPlus = (1/delta)*simp(F3.',delta_s);

            F4 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleLowerMidPoint^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleLowerMidPoint^2 + wireRadius^2);
            psi_nPlus_mMinus = (1/delta)*simp(F4.',delta_s);

            F5 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleLowerMidPoint^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleLowerMidPoint^2 + wireRadius^2);
            psi_nMinus_mMinus = (1/delta)*simp(F5.',delta_s);

            matrix = factor2*(psi_nPlus_mPlus-psi_nMinus_mPlus-psi_nPlus_mMinus+psi_nMinus_mMinus);
            esum = esum + matrix*I(n-1);
        end
        efield_Y(q,p) = abs(-esum/0.001/1.414);
        angle_efield_Y1(q,p) = rad2deg(angle(-esum));
    end
end

for q = 1:length(rho_vector)
    for p = 1:length(z_vector2)
        basisPoint_z(noOfBasisPoint+1) = z_vector2(p);
        basisPoint_z(noOfBasisPoint+2) = z_vector2(p);
        basisPoint_z(noOfBasisPoint+3) = z_vector2(p);
        testingDipoleLowerPoint_x = -testingDipoleHalfLength;
        tesingDipoleUpperPoint_x = testingDipoleHalfLength;
        testingDipoleMainPoint_x = 0;
        testingDipoleLowerMidPoint = (testingDipoleLowerPoint_x+testingDipoleMainPoint_x)/2;
        testingDipoleUpperMidPoint = (testingDipoleMainPoint_x+tesingDipoleUpperPoint_x)/2;
        
        esum = 0;
        for n = 2:noOfBasisPoint-1
            F2 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleUpperMidPoint^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleUpperMidPoint^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nPlus_mPlus = (1/delta)*simp(F2.',delta_s);

            F3 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleUpperMidPoint^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleUpperMidPoint^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nMinus_mPlus = (1/delta)*simp(F3.',delta_s);

            F4 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleLowerMidPoint^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleLowerMidPoint^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nPlus_mMinus = (1/delta)*simp(F4.',delta_s);

            F5 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleLowerMidPoint^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleLowerMidPoint^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nMinus_mMinus = (1/delta)*simp(F5.',delta_s);

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
plot(z_vector2,efield_rho(1,:),'r','LineWidth',2);hold on;
plot(z_vector2,efield_rho(2,:),'b','LineWidth',2);hold on;
plot(z_vector2,efield_rho(3,:),'m','LineWidth',2);hold on;


xlabel('z/\lambda','Fontsize',14);
ylabel('E_\rho in V/m','Fontsize',14);
legend('\rho/\lambda=0.01','\rho/\lambda=0.02','\rho/\lambda=0.05','Location','north','Orientation','vertical','Fontsize',14);


z_vector2 = z_vector2(2:length(z_vector2));
fig5 = figure('Color','w');
plot(z_vector2,angle_efield_Y1(1,2:length(angle_efield_Y1)),'r','LineWidth',2);hold on;
plot(z_vector2,angle_efield_Y1(2,2:length(angle_efield_Y1)),'b','LineWidth',2);hold on;
plot(z_vector2,angle_efield_Y1(3,2:length(angle_efield_Y1)),'m','LineWidth',2);


xlabel('z/\lambda','Fontsize',14);
ylabel('Phase of E_\rho (in degrees)','Fontsize',14);
ylim([-180 180]);
legend('\rho/\lambda=0.01','\rho/\lambda=0.02','\rho/\lambda=0.05','Location','north','Orientation','vertical','Fontsize',14);
toc;


function q=simp(y,dx)
N=length(y);

mul=ones(1,N);
mul(2:2:N-1)=4;
mul(3:2:N-2)=2;

q=(dx/3)*sum(y.'.*mul);
end