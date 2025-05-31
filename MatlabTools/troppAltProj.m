function F=troppAltProj(d,n,mu)
%troppAltProj
% Code to numerically approximate optimal line packings.
%
% Implementation of the algorithm in "Designing Structured Tight Frames 
% Via an Alternating Projection Method" by Joel A. Tropp, 
% Inderjit S. Dhillon, Robert W. Heath, Jr., and Thomas Strohmer and
% Tropp's dissertation "Topics in Sparse Approximation."
%
% INPUTS:
% d = dimension of vectors
% N = number of vectors
% mu = desired coherence to hit
%
% OUTPUTS:
% F = synthesis matrix of line packing
%
% For more information, see https://www.math.colostate.edu/~king/GameofSloanes.html
% and the paper "Game of Sloanes: Best known packings in complex projective space"
% by John Jasper, Emily J. King, and Dustin G. Mixon

% Initial Matrix (Algorithm 7.4 in Tropp's Thesis)
tau=0.9;
T=10000;
done=0;
while done == 0
    X=randn(d,1)+j*randn(d,1);
    X=X/norm(X);
    t=0;
    while t<T
        x=randn(d,1)+j*randn(d,1);
        x=x/norm(x);
        cor=max(max(abs(X'*x)));
        if cor<tau
            X(:,end+1)=x;
            t=0;
            if size(X,2)==n
                t=T;
            end
        else
            t=t+1;
        end
    end
    if size(X,2)==n
        done=1;
    end
end


% Alternating Projections (Algorithm 7.1 in Tropp's theis)
G=X'*X;
T=30000;
for t=1:T
    % project onto H(mu) by Proposition 7.6 in Tropp's theis
    for ii=1:n
        for jj=1:n
            if ii==jj
                G(ii,jj)=1;
            elseif abs(G(ii,jj))>mu
                G(ii,jj)=G(ii,jj)/abs(G(ii,jj))*mu;
            end
        end
    end
    G=(G+G')/2;
    % project onto G by Proposition 7.3 in Tropp's theis
    [V D] = eig(G);
    [dd,ind] = sort(diag(D),'descend');
    Vs = V(:,ind);
    first=dd(1:d);
    if sum(first)>n
        while sum(first)>n
            k=sum(first>0);
            gammapart=(sum(first)-n)/k;
            first=(first-gammapart).*((first-gammapart)>0);
        end
    else
        first=first+(n-sum(first))/d;
    end
    newD=zeros(size(D));
    newD(1:d,1:d)=diag(first);
    G=Vs*newD*Vs';
end

% finalize
D=diag(diag(G));
G=inv(D^(1/2))*G*inv(D^(1/2));
G=(G+G')/2;
[V D] = eig(G);
[dd,ind] = sort(diag(D),'descend');
Ds = D(ind,ind);
Vs = V(:,ind);
sqrtDs=sqrt(Ds);
Dnew=zeros(d,n);
Dnew(1:d,1:d)=sqrtDs(1:d,1:d);
F=Dnew*Vs';
for ii=1:n
    F(:,ii)=F(:,ii)/norm(F(:,ii));
end