%RWG3 FREQUENCY LOOP 
%   Calculates the impedance matrix using function IMPMET
%   and solves MoM equations
%   Uses the mesh file from RWG2, mesh2.mat, as an input.
%   Includes three additional impedance matrixes:

%   Dielectric-to-dielectric impedance matrix       ZDD
%   Dielectric-to-metal impedance matrix            ZDS 
%   Metal-to-dielectric impedance matrix            ZSD 
% 
%   The following parameters need to be specified prior to 
%   calculations:
%   
%   Number of frequency steps       NumberOfSteps
%   Lower frequency                 FreqStart
%   Upper frequency                 FreqStop
%   Dielectric constant (SI)        epsilon_
%   Magnetic permeability (SI)      mu_
%   Relative dielectric constant of
%   the dielectric material         epsilon_R
%
%   The algorithm of Chapter 10 uses the transpose impedance matrix
%   Z.' instead of Z for the metal structure. To better see what is 
%   the physical difference between Z and Z.' one should again return
%   to the basic paper Rao, Wilton, Glisson, IEEE Trans. Antennas and 
%   Propagation, vol. AP-30, No 3, pp. 409-418, 1982. It can be shown 
%   that Z utilizes 1 integration point for the outer integral (m) and 
%   9 integration points for the inner integral (n) whereas its transpose 
%   implies 9 integration points for the outer integral and 1 point for 
%   the inner integral. The difference between Z and Z.' is really 
%   insignificant for the pure metallic antennas. Either Z or Z.' can 
%   be used.    
%   However, for the metal-dielectric structures, Z.' may be preferred 
%   since it better matches other mixed integrals (metal-dielectric).  
%
%   Copyright 2002 AEMM. Revision 2002/03/23 
%   Chapter 10

clear all
%Load the data
load('mesh2');

%Frequency series parameters
NumberOfSteps=21;
FreqStart   =1.0e9;      %in Hz
FreqStop    =6.0e9;      %in Hz
step=(FreqStop-FreqStart)/(NumberOfSteps-1);

%EM parameters
epsilon_    =8.854e-012;
epsilon_R   =1.0;
mu_         =1.257e-006;
%speed of light 
c_=1/sqrt(epsilon_*mu_);
%free-space impedance 
eta_=sqrt(mu_/epsilon_);

%Contemporary variables - metal impedance matrix
for m=1:EdgesTotal
    RHO_P(:,:,m)=repmat(RHO_Plus(:,m),[1 9]);   %[3 9 EdgesTotal]
    RHO_M(:,:,m)=repmat(RHO_Minus(:,m),[1 9]);  %[3 9 EdgesTotal]
end

%Contemporary variables - dielectric/metal impedance matrix
DP      =find(t(4,:)==0);
M       =length(DP);        %first M triangles form the dielectric plate
delta   =[0;0;1e-12];       %introduce to avoid singularity 0/0

Middle=[0; 0; h/2];
for m=1:M
    N=t(1:3,m);
    Point(:,m)=Center(:,m)      +Middle; 
    IMT(:,:,m)=Center_(:,:,m)   +repmat(Middle,[1,9]);
end

%Frequency series
for FF=1:NumberOfSteps    
    tic; FF
    f(FF)       =FreqStart+step*(FF-1);
    omega       =2*pi*f(FF);
    k           =omega/c_;
    K           =j*k;
  
    Constant1   =mu_/(4*pi);
    Constant2   =1/(j*4*pi*omega*epsilon_);
    Factor      =1/9;    
    FactorA     =Factor*(j*omega*EdgeLength/4)*Constant1;
    FactorFi    =Factor*EdgeLength*Constant2;
    FactorA     =FactorA.';
    FactorFi    =FactorFi.';
    
    %Metal-to-metal impedance matrix (SS)
    ZSS=  impmet( EdgesTotal,TrianglesTotal,...
            EdgeLength,K,...
            Center,Center_,...
            TrianglePlus,TriangleMinus,...
            RHO_P,RHO_M,...
            RHO__Plus,RHO__Minus,...
            FactorA,FactorFi);   
    ZSS=ZSS.'; 
    toc; tic;
    
    %Dielectric-to-dielectric impedance matrix (DD) 
    ZDD = zeros(M,M)+j*zeros(M,M);
    for m=1:M
        OP      =Point(:,m)+delta;
        IP      =Point;
        E       =point_(OP,K,k,Constant2,Area(1:M),h,IP);
        ZDD(m,:)=E(3,:);        %Non-diagonal terms
        IP      =IMT(:,:,m);
        E       =point_(OP,K,k,Constant2,ones(1,9)*Area(m)/9,h,IP);
        ZDD(m,m)=sum(E(3,:));   %Diagonal terms
    end    
    ZDD=j*omega*epsilon_*(epsilon_R-1)*ZDD; %Non-symmetric matrix since areas are different
    toc; tic; 
     
    %Dielectric-to-metal impedance matrix (DS) 
    ZDS = zeros(EdgesTotal,M)+j*zeros(EdgesTotal,M);
    for m=1:EdgesTotal
        OPPlus  =Center(:,TrianglePlus(m))+delta;
        OPMinus =Center(:,TriangleMinus(m))+delta;
        EP      =point_(OPPlus,K,k,Constant2,Area(1:M),h,Point);
        EM      =point_(OPMinus,K,k,Constant2,Area(1:M),h,Point);
        ScalarPlus  =sum(EP.* repmat(RHO_Plus(:,m), [1 M]));
        ScalarMinus =sum(EM.* repmat(RHO_Minus(:,m),[1 M]));
        ZDS(m,:)    =EdgeLength(m)*(ScalarPlus/2+ScalarMinus/2);    
    end
    toc; tic;
       
    %Metal-to-dielectric impedance matrix (SD) 
    ZSD = zeros(M,EdgesTotal)+j*zeros(M,EdgesTotal);
    C1=1/(2*epsilon_)/(-j*omega);
    C2=j*omega*epsilon_*(epsilon_R-1);
    C=C1*C2;
    for m=1:M         
        Q=sum(abs([Center(1,m)-Center(1,:); Center(2,m)-Center(2,:)]));
        T=find(Q<1.e-9); %Finds only bottom or bottom and top triangles       
        for q=1:length(T)
            Plus    =find(TrianglePlus-T(q)==0);
            Minus   =find(TriangleMinus-T(q)==0);
            Ind     =(-1)^(q+1);
            for k=1:length(Plus)
                n=Plus(k);
                Charge=EdgeLength(n)/Area(T(q));
                ZSD(m,n)=ZSD(m,n)+C*Ind*Charge;
            end
            for k=1:length(Minus)
                n=Minus(k);
                Charge=-EdgeLength(n)/Area(T(q));
                ZSD(m,n)=ZSD(m,n)+C*Ind*Charge; 
            end 
        end
    end
    toc; tic;
  
    for m=1:EdgesTotal    
        V(m)=0;
    end
    
    %This is the bottom probe feed
    Index=find(EdgeIndicator==1);
    
    %This is the center probe feed 
    %(use for the dielectric example of Section 10.9)
    %for m=1:EdgesTotal
    %    ECenter(:,m)=1/2*(p(:,Edge_(1,m))+p(:,Edge_(2,m)));
    %end
    %[Y,INDEX]=sort(abs(ECenter(3,:)-repmat(-h/2,[1,EdgesTotal])));
    %Index=INDEX(1);

    %Solution of MoM equations
    V(Index)=1*EdgeLength(Index);    
    %Total impedance matrix
    Z=[ZSS -ZDS; ZSD ZDD-eye(M,M)];
    %Total voltage vector
    V(EdgesTotal+1:EdgesTotal+M)=0;
    
    I=Z\V.';
    toc; tic;
    
    CURRENT(:,FF)=I(:);
    %Impedance
    GapCurrent(FF)  =sum(I(Index).*EdgeLength(Index)');
    GapVoltage(FF)  =mean(V(Index)./EdgeLength(Index));
    Impedance(FF)   =GapVoltage(FF)/GapCurrent(FF);
    FeedPower(FF)   =1/2*real(GapCurrent(FF)*conj(GapVoltage(FF)));    
    Imp             =Impedance(FF)
end

%Save result
FileName='current.mat'; 
save(FileName, 'f','NumberOfSteps','FreqStart','FreqStop','step',...
                'omega','mu_','epsilon_','c_', 'eta_',...
                'CURRENT','GapCurrent','GapVoltage','Impedance','FeedPower','M','h','Point','Index');            