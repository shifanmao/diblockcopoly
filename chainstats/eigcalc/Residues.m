function out=Residues(K,R,ORDEig,ORDL,ResLayer,ORDMU)
%function [resi resia]=gl1l2m(R,k,d,L1,L2,M,ORDL,ResLayer)
% use asymptotic residual to figure out whether calculate the residual
% using continued fraction or stay with the asymptotic form for small k
% limit

d=3;
NR=length(R);           % number of eigenvalues (roots)
ResThreshold=1e-12;     % threshold to go from small k asymptot to matrix method
ImagThreshold=1e-8;     % threshold to go from small k asymptot to matrix method

cutoff=1e-11; % chosen by looking at graph
out=zeros(ORDL,ORDL,ORDMU,ORDEig)*NaN;

for lam1=0:(ORDL-1)
    for lam2=0:(ORDL-1)
        for mu=0:min([lam1,lam2,ORDMU-1])
	    smallK=SmallAsympRes(K,NR,d,lam1,lam2,mu);            
	    out(lam1+1,lam2+1,mu+1,:)=smallK;

            % apply cutoff and evalute using continued fraction
            cutoffindex = find(abs(smallK)>cutoff);
	    out(lam1+1,lam2+1,mu+1,cutoffindex)=...
	      gl1l2m(K,R(cutoffindex),d,lam1,lam2,mu,ResLayer);
        end
    end
end

end

function resi=gl1l2m(k,R,d,L1,L2,M,ResLayer)
% get the residues for all roots given in R using recursive relation for
% derivative

NR=length(R);           % number of eigenvalues (roots)	  
resi=zeros(NR,1);
for iEig=1:1:NR
		  
    L=sort([L1 L2]);
    L1=L(1);
    L2=L(2);

    AL=zeros(ResLayer,1);
    WP=zeros(ResLayer,1);
    WM=zeros(ResLayer,1);
    dJp=zeros(ResLayer,1);
    dJm=zeros(ResLayer,1);
    
    n=ResLayer-1;
    AL(ResLayer)=sqrt((n-M)*(n+M+d-3)/(2*n+d-2)/(2*n+d-4));
    WP(ResLayer)=1/(R(iEig)+n*(n+d-2));
    dJp(ResLayer)=1;

    for n=ResLayer-2:-1:M    
        PL=R(iEig)+n*(n+d-2);        
        AL(n+1)=sqrt((n-M)*(n+M+d-3)/(2*n+d-2)/(2*n+d-4));
        
        WP(n+1)=1/(PL+WP(n+2)*(AL(n+2)*k)^2);
        dJp(n+1)=1-dJp(n+2)*(AL(n+2)*k*WP(n+2))^2;
    end

    n=M;
    PL=R(iEig)+n*(n+d-2);            
    WM(n+1)=1/PL;
    dJm(n+1)=1;
    
    for n=(M+1):1:L1
        PL=R(iEig)+n*(n+d-2);    
        AL(n+1)=sqrt((n-M)*(n+M+d-3)/(2*n+d-2)/(2*n+d-4));
        
        WM(n+1)=1/(PL+WM(n)*(AL(n+1)*k)^2);
        dJm(n+1)=1-dJm(n)*(AL(n+1)*k*WM(n))^2;                
    end

    if L1 > M
        wl1m=1/(-(AL(L1+1)*k*WM(L1))^2*dJm(L1)+1-(AL(L1+2)*k*WP(L1+2))^2*dJp(L1+2));
    else
        wl1m=1/dJp(L1+1);        
    end
    
    if L1 == L2
        resi(iEig)=wl1m;
    else
        resi(iEig)=wl1m*prod(1i*k*AL((L1+1+1):(L2+1)).*WP((L1+1+1):(L2+1)));
    end
    
end
    
end

function resia=SmallAsympRes(K,NR,d,L1,L2,M)

% calculate the residue using small k asymptot
resia=zeros(NR,1);

for L=0:1:(NR-1)
    
    if L==L1
        CL1=1;
    elseif L > L1
        CL1=1;
        for n=L1:(L-1)
            AL=sqrt((n+1-M)*(n+1+M+d-3)/(2*n+d)/(2*n+d-2));
            CL1=CL1*AL/(L*(L+d-2)-n*(n+d-2));
        end
        CL1=CL1*(-1i*K)^(L-L1);        
    elseif L < L1
        CL1=1;
        for n=(L+1):L1
            AL=sqrt((n-M)*(n+M+d-3)/(2*n+d-2)/(2*n+d-4));
            CL1=CL1*AL/(L*(L+d-2)-n*(n+d-2));
        end
        CL1=CL1*(-1i*K)^(L1-L);                
    end
    
    if L==L2
        CL2=1;
    elseif L > L2
        CL2=1;
        for n=L2:(L-1)
            AL=sqrt((n+1-M)*(n+1+M+d-3)/(2*n+d)/(2*n+d-2));
            CL2=CL2*AL/(L*(L+d-2)-n*(n+d-2));
        end
        CL2=CL2*(-1i*K)^(L-L2);        
    elseif L < L2
        CL2=1;
        for n=(L+1):L2
            AL=sqrt((n-M)*(n+M+d-3)/(2*n+d-2)/(2*n+d-4));
            CL2=CL2*AL/(L*(L+d-2)-n*(n+d-2));
        end
        CL2=CL2*(-1i*K)^(L2-L);                
    end
    
    resia(L+1)=CL1*CL2;
    
end

end
