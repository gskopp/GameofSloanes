function [ Phi, mu ] = GrassPack( d, N, coh, F )  %remove last input for script
%GrassPack 
% Code to numerically approximate optimal line packings.
% Modification of the code from 
% http://www.ece.mcmaster.ca/~davidson/pubs/Flexible_codebook_design.html,
% which is the implementation of "Flexible codebook design for limited 
% feedback systems via sequential smooth optimization on the Grassmannian
% manifold" by Ahmed Medra and Timothy N. Davidson
%
% INPUTS:
% d = dimension of vectors
% N = number of vectors
% coh = desired coherence to beat
%
% OUTPUTS:
% Phi = synthesis matrix of line packing
% mu = coherence of Phi
%
% Created on July 11, 2018
%
% For more information, see https://www.math.colostate.edu/~king/GameofSloanes.html
% and the paper "Game of Sloanes: Best known packings in complex projective space"
% by John Jasper, Emily J. King, and Dustin G. Mixon



% clc
% close all
% clear all
global M
global Mt
global Qpre
global skip
global power
global WM


Mt=d;
M=1; % doing line packings
skip=0;
Qpre=N;
QLeft=Qpre;

%x=[genorth(Mt,M,QLeft)]; %uncomment for script
x= F; % load another packing
WM=genweight(Mt,M,Qpre);
tol=10^-4;
mag=1;
ma=10;
power=2;

for jj=1:2000%200
    mag=mag*10;
    [x fval]=orthcon(x,@objectivefn,mag);
    [x fval]=orthcon(x,@objectivefn2,ma);
    if fval<tol
        break
    end
    power=power+2;
    if power==20
        power=2;
        mag=1;
    end
    if mod(jj,5)==0
        %close all
        mu=DispRes(x);
        Phi=x;
        if mu<coh
            break
        end
    end
%     if mod(jj,5)==0
%         close all
%         ss=input('Do you wish to continue? Press 1 to stop or Enter to continue');
%         if ss==1
%             break
%         end
%     end
end

%DispRes(x)

%==========================================================================
%==========================================================================
function [x fval]=orthcon(x,objectivefn,mag)
%====init=====
global  M Mt
global Qpre
global skip
maxiter=20;
xbest=x;
fbest=objectivefn(x);
direction=[zeros(1,skip) ones(1,(Qpre-skip))];
rate=0;
%=========================================================================
for kkk=1:maxiter
gamma=ones(1,Qpre);
Dx=mag*(NumGradI(objectivefn,(x)));
Z=[];
ttrace=[];
cont=[];
  
for jj=1:Qpre
    f=x(:,(jj-1)*M+1:jj*M);
    h=f*f';
    X=Dx(:,(jj-1)*M+1:jj*M);
    I1=eye(size(h));
    con1=(h-I1);
    z1=-(I1-h)*(X);
    cont=[cont con1];
    Z=[Z z1];
    ttrace=[ttrace trace(z1'*z1)];
end

XX=x+Z;
xnn=x;
xn=xnn;
dd=feval(objectivefn,x);
for jj=(skip+1):Qpre
    xn=xnn;
    XX1=XX(:,(jj-1)*M+1:jj*M);

    [Q1 R1]=qr(XX1);
    xn(:,(jj-1)*M+1:jj*M)=Q1(:,1:M);
    QRfac=feval(objectivefn,xn);
    ii=0;

    QRfacinit=QRfac;
    fact=1;
    while (dd-QRfac)>= gamma(jj)*ttrace(jj)*0 && ii<50
    gamma(jj)=2*gamma(jj);
    temp=x(:,(jj-1)*M+1:jj*M)+gamma(jj)*Z(:,(jj-1)*M+1:jj*M);
    [Q1 R1]=qr(temp);
    xn(:,(jj-1)*M+1:jj*M)=Q1(:,1:M);
    QRfac=feval(objectivefn,xn);
    if QRfac>QRfacinit
        break
    end
    QRfacinit=QRfac;
    ii=ii+1;
    fact=2;
    end
    gamma(jj)=gamma(jj)/fact;
    temp=x(:,(jj-1)*M+1:jj*M)+gamma(jj)*Z(:,(jj-1)*M+1:jj*M);
    [Q1 R1]=qr(temp);
    xn(:,(jj-1)*M+1:jj*M)=Q1(:,1:M);
    QRfac=feval(objectivefn,xn);
    ii=0;
    fact=1;
    while (dd-QRfac)<0.5*gamma(jj)*ttrace(jj)*0 && ii<30
    gamma(jj)=gamma(jj)*0.25;
    temp=x(:,(jj-1)*M+1:jj*M)+gamma(jj)*Z(:,(jj-1)*M+1:jj*M);
    [Q1 R1]=qr(temp);
    xn(:,(jj-1)*M+1:jj*M)=Q1(:,1:M);
    QRfac=feval(objectivefn,xn);
    ii=ii+1;

    end
    gamma(jj)=gamma(jj)/fact;
    temp=x(:,(jj-1)*M+1:jj*M)+gamma(jj)*Z(:,(jj-1)*M+1:jj*M);
    [Q1 R1]=qr(temp);
    xnn(:,(jj-1)*M+1:jj*M)=Q1(:,1:M);
    dd=feval(objectivefn,xnn);
end
xn=x;
for jj=(skip+1):Qpre
        XX1=xnn(:,(jj-1)*M+1:jj*M);
        [Q1 R1]=qr(XX1);
        xn(:,(jj-1)*M+1:jj*M)=Q1(:,1:M);
end
x=xn;
fval=feval(objectivefn,x);
disp(['iteration number = ',num2str(kkk),'  Function value = ',num2str(fval)]);

%=========================================================================
mm=fval;
if mm<fbest
%     disp('fbest changed');
    fbest=mm;
    xbest=x;
    rate=0;

elseif mm>=fbest
%     disp(['fbest unchanged',num2str(fbest)]);
%     mag=max(mag*0.95,0.5)
    rate=rate+1;
    if rate==2
        x=xbest;
       mag=mag/10;
    end
    if rate==4
        x=xbest;
       mag=mag*100;
    end
    if rate==6
        x=xbest;
       mag=mag/1000;
    end

    if rate==10
        x=xbest;
        break;
    end
end
%==========================================================================    

nn=isnan(x);
if nn>0 
    disp('Nan')
%     [x fval]=OrthogonalAlgGen3()
    return
end
end
disp('Maximum iteration reached!')
fval=feval(objectivefn,xbest);
DispRes(xbest)
x=xbest;
%=========================================================================
beep
disp('====================================================================')
return    
   


%=========================================================================
function [Dx] = NumGradI(fcn,x)
global M
% this code is modified from another one.
%-------------------------------------------------------------
f0 = feval(fcn,x);
t=zeros(size(x));
delta = 1e-4;
gradR=t;
gradI=t;
for ii=1:1:size(x,1)
    for kk=1:1:size(x,2)
        z=t;z(ii,kk)=1*delta;
        xx=x+z;
        f1=feval(fcn,xx);
        gradR(ii,kk)=(f1-f0)/delta;
    end
end
for ii=1:1:size(x,1)
    for kk=1:1:size(x,2)
        z=t;z(ii,kk)=1i*delta;
        xx=x+z;
        f1=feval(fcn,xx);
        gradI(ii,kk)=(f1-f0)/delta;
    end
end
Dx=gradR+1i*gradI;

    

%==========================================================================
function f=objectivefn(F)
global M
global Qpre
global power
global WM
Mt=size(F,1);
Qpred=Qpre;
cor=1-(Mt-1)*Qpred/Mt/(Qpred-1);
if Qpre>Mt^2
    cor=0;
end
fac=sqrt(cor);
faclog=((1+det(fac*fac))*log(1+det(fac*fac)));
% ee2=[];
MM=F'*F;
MM=(1+(conj(MM).*MM)).*(log(1+conj(MM).*MM));
MM=(MM-faclog).*WM;
ee2=reshape(MM,numel(MM),1);
% for ii=1:Qpred
%     for jj=(skip+1):Qpred
%         if ii==jj
%             continue
%         end
%         f1=F(:,(ii-1)*M+1:ii*M);
%         f2=F(:,(jj-1)*M+1:jj*M);
%         d2=f1'*f2;
%         ee2=[ee2 ((1+det(d2'*d2))*log(1+det(d2'*d2)))-faclog];
%     end
% end
f=(mean((ee2).^power));
%=========================================================================
function cohmx=DispRes(F)
% global M

Qpre=size(F,2);

cohmx=max(max(abs(F'*F)-eye(Qpre)));
disp(['Coherence = ',num2str(cohmx)])

% ddd=[];
% Qpre=size(F,2);
% 
% for ii=1:Qpre
%     for jj=ii:Qpre
%         if ii==jj
%             continue
%         end
%         f1=F(:,(ii-1)*M+1:ii*M);
%         f2=F(:,(jj-1)*M+1:jj*M);
%         d=abs(f2'*f1);
%         %d=sqrt(1-d^2);
%         ddd=[ddd d];
%      
%     end
% end
% % chmin=min(ddd);
% % disp(['Minimum chordal distance = ',num2str(chmin)])
% cohmx=max(ddd);
% disp(['Coherence = ',num2str(cohmx)])
% % figure
% % hist((ddd))
% % xlabel('Distance Values')
% % ylabel('Distance Distribution')

%==========================================================================
function q=genorth(Mt,M,Qpre)
q=[];
for ii=1:Qpre
    q=[q orth(randn(Mt,M)+1i*randn(Mt,M))];
end
%==========================================================================
%-------------------------------------------------------------------------
function f=objectivefn2(F)
global M
global Qpre
global WM
Mt=size(F,1);
Qpred=Qpre;
cor=1-(Mt-1)*Qpred/Mt/(Qpred-1);
if Qpre>Mt^2
    cor=0;
end
fac=sqrt(cor);
faclog=((1+det(fac*fac))*log(1+det(fac*fac)));
MM=F'*F;
MM=(1+(conj(MM).*MM)).*(log(1+conj(MM).*MM));
MM=(MM-faclog).*WM;
ee=reshape(MM,numel(MM),1);
% for ii=1:Qpre
%     for jj=(skip+1):Qpre
%         if ii==jj
%             continue
%         end
%         f1=F(:,(ii-1)*M+1:ii*M);
%         f2=F(:,(jj-1)*M+1:jj*M);
%         d=f2'*f1;
%         ee=[ee ((1+det(d'*d))*log(1+det(d'*d)))-faclog];
%     end
% end
f=norm(ee,'inf');
%==========================================================================
%==========================================================================
function correq=bestknown(dim,N)
cor=(dim-1)*N/dim/(N-1);
correq=sqrt(cor);
%==========================================================================
function WM=genweight(Mt,M,Qpre)
WM=zeros(Mt,M*Qpre);
for ii=1:Qpre
    for jj=ii:Qpre
       WM(ii,jj)=1;
    end
end
WM=1-WM;


