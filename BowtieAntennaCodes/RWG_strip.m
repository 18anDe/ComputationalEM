clear;close;clc;

tic;

%% Antenna Parameters

% Change the frequency of the antenna and the antenna structure for
% different lambda values

lambda = 3;       %Change this for different lambda values
mult = 0.75;           %Change accordingly for bowtie antenna
eps0 = 8.8542e-12;
mu0 = pi*4e-7;
c0 = 1/sqrt(mu0*eps0);
eta0 = sqrt(mu0/eps0);
f0 = c0/lambda;
omega = 2*pi*f0;
k = omega/c0;
K = 1j*k;
noOfSegments = 100;        % To match with the strip current the no of segments are taken high
stripLength = lambda/2;
%stripWidth = (0.005*lambda)/0.25;      % DipoleAntenna
stripWidth = lambda/40;

%% Antenna Structure

load('bowtie.mat');         % Change this for different lambda values
p(3,:) = 0;                 % to convert 2D to 3D
TrianglesTotal=length(t);

%% Area and Center of Triangles

Area = zeros;
Center = zeros(3,TrianglesTotal);
for m=1:TrianglesTotal
   N=t(1:3,m);
   Vec1=p(:,N(1))-p(:,N(2));
   Vec2=p(:,N(3))-p(:,N(2));
   Area(m) =norm(cross(Vec1,Vec2))/2;
   Center(:,m)=1/3*sum(p(:,N),2);
end

%% Common Edges for adjacent triangles

Edge_=[];
TrianglePlus = zeros;
TriangleMinus = zeros;
n=0;
for m=1:TrianglesTotal
    N=t(1:3,m);
    for r=m+1:TrianglesTotal
        M=t(1:3,r);      
        a=1-all([N-M(1) N-M(2) N-M(3)]);
        if(sum(a)==2) %triangles m and k have common edge
            n=n+1;
            Edge_=[Edge_ M(find(a))]; 
            TrianglePlus(n)=m;
            TriangleMinus(n)=r; 
        end
    end
end
EdgesTotal=length(Edge_);

%% Edge Length

EdgeLength = zeros;
for m=1:EdgesTotal
   EdgeLength(m)=norm(p(:,Edge_(1,m))-p(:,Edge_(2,m)));
end

%% Midpoints of nine subtriangles

Center_ = zeros(3,9,TrianglesTotal);
for m=1:TrianglesTotal
    n1=t(1,m);
    n2=t(2,m);
    n3=t(3,m); 
    MP=Center(:,m);
    r1=    p(:,n1);
    r2=    p(:,n2);
    r3=    p(:,n3);
    r12=r2-r1;
    r23=r3-r2;
    r13=r3-r1;
    C1=r1+(1/3)*r12;
    C2=r1+(2/3)*r12;
    C3=r2+(1/3)*r23;
    C4=r2+(2/3)*r23;
    C5=r1+(1/3)*r13;
    C6=r1+(2/3)*r13;
    a1=1/3*(C1+C5+r1);
    a2=1/3*(C1+C2+MP);
    a3=1/3*(C2+C3+r2);
    a4=1/3*(C2+C3+MP);
    a5=1/3*(C3+C4+MP);
    a6=1/3*(C1+C5+MP);
    a7=1/3*(C5+C6+MP);
    a8=1/3*(C4+C6+MP);
    a9=1/3*(C4+C6+r3);
    Center_(:,:,m)=...
        [a1 a2 a3 a4 a5 a6 a7 a8 a9];
end

%% RHO vectors for plus and minus triangles and the subtriangles

%PLUS

RHO_Plus = zeros(3,EdgesTotal);
RHO__Plus = zeros(3,9,EdgesTotal);

for m=1:EdgesTotal
    NoPlus=TrianglePlus(m);
    n1=t(1,NoPlus);
    n2=t(2,NoPlus);
    n3=t(3,NoPlus); 
    if n1~=Edge_(1,m) && n1~=Edge_(2,m)
        NODE=n1; 
    end
    if n2~=Edge_(1,m) && n2~=Edge_(2,m)
        NODE=n2; 
    end
    if n3~=Edge_(1,m) && n3~=Edge_(2,m) 
        NODE=n3; 
    end
    FreeVertex=p(:,NODE);
    
    RHO_Plus(:,m)   =+Center(:,NoPlus)-FreeVertex;
    %Nine rho's of the "plus" triangle
    RHO__Plus(:,:,m)  =...
        +Center_(:,:,NoPlus)-repmat(FreeVertex,[1 9]);
end

%MINUS

RHO_Minus = zeros(3,EdgesTotal);
RHO__Minus = zeros(3,9,EdgesTotal);

for m=1:EdgesTotal
    NoMinus=TriangleMinus(m);
    n1=t(1,NoMinus);
    n2=t(2,NoMinus);
    n3=t(3,NoMinus); 
    if n1~=Edge_(1,m) && n1~=Edge_(2,m)
        NODE=n1; 
    end
    if n2~=Edge_(1,m) && n2~=Edge_(2,m)
        NODE=n2; 
    end
    if n3~=Edge_(1,m) && n3~=Edge_(2,m) 
        NODE=n3; 
    end
    FreeVertex=p(:,NODE);
    
    RHO_Minus(:,m)   =-Center(:,NoMinus) +FreeVertex;
    %Nine rho's of the "minus" triangle
    RHO__Minus(:,:,m)=...
        -Center_(:,:,NoMinus)+repmat(FreeVertex,[1 9]);
end

%% Constants

Constant1 = mu0/(4*pi);
Constant2 = 1/(1j*4*pi*omega*eps0);
Factor = 1/9;
FactorA = Factor*(1j*omega*EdgeLength.'/4)*Constant1;
FactorFi = Factor*EdgeLength.'*Constant2;

%% RHO vectors

RHO_P = zeros(3,9,EdgesTotal);
RHO_M = zeros(3,9,EdgesTotal);
for m=1:EdgesTotal
    RHO_P(:,:,m)=repmat(RHO_Plus(:,m),[1 9]);   
    RHO_M(:,:,m)=repmat(RHO_Minus(:,m),[1 9]);  
end

%% Impedance Matrix

Z = impmet( EdgesTotal,TrianglesTotal,EdgeLength,K,Center,Center_,...
            TrianglePlus,TriangleMinus,RHO_P,RHO_M,RHO__Plus,RHO__Minus,...
            FactorA,FactorFi);  

%% Feeding edge, Voltage vector and solving MoM equations

FeedPoint=[0; 0; 0];

V = zeros;
Distance = zeros(3,EdgesTotal);

for m=1:EdgesTotal
    V(m)=0;
    Distance(:,m)=0.5*sum(p(:,Edge_(:,m)),2)-FeedPoint;
end

[Y,INDEX]=sort(sum(Distance.*Distance));
Index=INDEX(1); 

V(Index)=1*EdgeLength(Index);    

I=Z\V.';

%% Antenna input Impedance

GapCurrent  =sum(I(Index).*EdgeLength(Index)');
GapVoltage  =mean(V(Index)./EdgeLength(Index));
Impedance   = GapVoltage/GapCurrent;
FeedPower   = 1/2*real(GapCurrent*conj(GapVoltage));

%% Current Density

CurrentNorm = zeros;
CurrentNorm_new = zeros;

for q=1:TrianglesTotal
    i=[0 0 0]';
    for m=1:EdgesTotal
        IE=I(m)*EdgeLength(m);
        if(TrianglePlus(m)==q)
            i=i+IE*RHO_Plus(:,m)/(2*Area(TrianglePlus(m)));
        end
        if(TriangleMinus(m)==q)
            i=i+IE*RHO_Minus(:,m)/(2*Area(TriangleMinus(m)));
        end
    end
    CurrentNorm(q)=abs(norm(i));
    CurrentNorm_new(1:2,q)=(i(1:2));

end

Jmax=max(CurrentNorm);
MaxCurrent=strcat(num2str(Jmax),'[A/m]')
%CurrentNorm1=CurrentNorm/max(CurrentNorm);

X1 = zeros;
Y1 = zeros;
Z1 = zeros;

for m=1:TrianglesTotal
    N=t(1:3,m);
    X1(1:3,m)=[p(1,N)]';
    Y1(1:3,m)=[p(2,N)]';
    Z1(1:3,m)=[p(3,N)]';      
end
C=repmat(CurrentNorm,3,1);

fig1 = figure();
h=fill3(X1, Y1, Z1, C); %linear scale
colormap jet(20);
brighten(0.5);
minValue = min(CurrentNorm);
maxValue = max(CurrentNorm);
clim([minValue, maxValue]);
colorbar;
h1 = colorbar;
ticks = linspace(minValue, maxValue, 20);  % Adjust minValue, maxValue, and numTicks based on your preferences
h1.Ticks = ticks;
%view(0, 90);
axis('equal');
rotate3d on

%% rwg 6

K1=20;
x0=min(p(1,:));
x1=max(p(1,:));
y0=min(p(2,:));
y1=max(p(2,:));

y = zeros;
X2 = zeros;
Y2 = zeros;

for n=1:K1+1
    y(n)=y0+(n-1)*(y1-y0)/K1;
    Dist=repmat([0 y(n) 0]',[1,TrianglesTotal])-Center;
    [dummy,Index]=min( sum(Dist.*Dist));
    X2(n)=CurrentNorm_new(1,Index); 
    Y2(n)=CurrentNorm_new(2,Index); 
end
yi=y0:(y1-y0)/noOfSegments:y1;
Xi = interp1(y,X2,yi,'cubic');
Yi = interp1(y,Y2,yi,'cubic');

%I_strip = Yi*stripWidth;    %Dipole Antenna
I_strip = Yi*stripWidth;     %Bowtie Antenna
I_strip(:,1) = []; I_strip(:,length(I_strip)) = [];
I_strip1 = I_strip*1e3;

fig2 = figure();
yyaxis left
plot(yi,[0 abs(I_strip1) 0],'b','LineWidth',1.4);
%ylim([0 12])
ylabel('Magnitude of Current (mA)','FontSize',14)
hold on;
yyaxis right
plot(yi,[rad2deg(angle(I_strip1(1))) rad2deg(angle(I_strip1)) rad2deg(angle(I_strip1(noOfSegments-1)))],'--','LineWidth',1.4);     %Phase plot not correct
ylim([-180 180])
ylabel('Phase of Current (Degree)','FontSize',14)
xlabel('Dipole length, m','FontSize',14)
xlim([-stripLength/2 stripLength/2])
grid on
legend('Magnitude','Phase','Location','northwest','FontSize',14);

%% Results data

Filename = 'mesh2.mat';
save(Filename, "Center","EdgeLength","EdgesTotal","TriangleMinus","TrianglePlus")

Filename = 'current.mat';
save(Filename,"stripWidth","noOfSegments","I_strip", "I_strip1",...
    "stripLength","lambda","k","omega","mu0","eps0","K","I","eta0")

toc;