clear;close;clc;

tic;
    
%% Antenna Parameters

% Change the frequency of the antenna and the antenna structure for
% different lambda values

lambda = 1.44;       %Change this for different lambda values
eps0 = 8.8542e-12;
mu0 = pi*4e-7;
c0 = 1/sqrt(mu0*eps0);
eta0 = sqrt(mu0/eps0);
f0 = c0/lambda;
omega = 2*pi*f0;
k = omega/c0;
K = 1j*k;
noOfSegments = 100;        % To match with the strip current the no of segments are taken high
antennaLength = 1;
%stripWidth = (0.005*lambda)/0.25;      % Dipole Antenna
antennaWidth = 0.01;                % Bowtie Antenna

%% Antenna Structure


load('mesh\bowtie.mat');         % Change this for different lambda values
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
noOfSteps = 200;
betaLStart = 5;
betaLStop = 290;
FreqStart = betaLStart*c0/(360*antennaLength);
FreqStop = betaLStop*c0/(360*antennaLength);
step = (FreqStop-FreqStart)/(noOfSteps-1);

%% RHO vectors
RHO_P = zeros(3,9,EdgesTotal);
RHO_M = zeros(3,9,EdgesTotal);
for m=1:EdgesTotal
    RHO_P(:,:,m)=repmat(RHO_Plus(:,m),[1 9]);   
    RHO_M(:,:,m)=repmat(RHO_Minus(:,m),[1 9]);  
end

FeedPoint=[0; 0; 0];
    V = zeros;Distance = zeros(3,EdgesTotal);
for FF = 1:noOfSteps
    FF
    f(FF) = FreqStart+step*(FF-1);
    omega = 2*pi*f(FF);
    k = omega/c0;
    K = 1j*k;
    L = 360*antennaLength*f/c0;
    
    Constant1 = mu0/(4*pi);
    Constant2 = 1/(1j*4*pi*omega*eps0);
    Factor = 1/9;
    FactorA = Factor*(1j*omega*EdgeLength.'/4)*Constant1;
    FactorFi = Factor*EdgeLength.'*Constant2;
    
    Z = impmet( EdgesTotal,TrianglesTotal,EdgeLength,K,Center,Center_,...
            TrianglePlus,TriangleMinus,RHO_P,RHO_M,RHO__Plus,RHO__Minus,...
            FactorA,FactorFi);  
      
    for m=1:EdgesTotal
        V(m)=0;
        Distance(:,m)=0.5*sum(p(:,Edge_(:,m)),2)-FeedPoint;
    end
    [Y,INDEX]=sort(sum(Distance.*Distance));
    Index=INDEX(1); 

    V(Index)=1*EdgeLength(Index);    

    I=Z\V.';
    
    CURRENT(:,FF) = I(:);
    GapCurrent(FF)  =sum(I(Index).*EdgeLength(Index)');
    GapVoltage(FF)  =mean(V(Index)./EdgeLength(Index));
    Impedance(FF)   = GapVoltage(FF)/GapCurrent(FF);
    FeedPower(FF)   = 1/2*real(GapCurrent(FF)*conj(GapVoltage(FF)));
    Imp = Impedance(FF);
end

a=figure;
plot(L, real(Impedance),'b','LineWidth',1.4);
%xlabel ('Frequency, GHz','FontSize',14)
xlabel ('Antenna Length, Degrees','FontSize',14)
ylabel('Input  resistance, Ohm','FontSize',14)
xlim([0 300])
%axis([0 8000e6 0 400])
%grid on
b=figure;
plot(L, imag(Impedance),'b','LineWidth',1.4);
%xlabel ('Frequency, GHz','FontSize',14)
xlabel ('Antenna Length, Degrees','FontSize',14)
ylabel('Input  reactance, Ohm','FontSize',14)
xlim([0 300])
%axis([0 8000e6 -250 150])
%grid on

ax = gca;
ax.FontSize = 14;           % Font size adjusted to 14