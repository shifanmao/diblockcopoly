function [smallK,largeK]=Residues(K,EigK,ORDEig,ORDL,NumLayer,ORDMU)
%function out=Residues(K,EigK,ORDEig,ORDL,NumLayer,ORDMU)
% To speed this function up you may not need to return all M values,
% you could separate out a ORDM
% Returns a ORDL x ORDL x ORDMU x ORDEig matrix
% Index (lam1+1, lam2+1, mu+1, l+1)
% I haved included the ORDMU separate from ORDL to save computation time

cutoff=10^(-11); % chosen by looking at graph

largeK=glm(K,EigK,ORDEig,ORDL,NumLayer,ORDMU);

smallK=zeros(ORDL,ORDL,ORDL,ORDEig)*NaN;
out=zeros(ORDL,ORDL,ORDMU,ORDEig)*NaN;

for lam1=0:ORDL-1
    for lam2=0:ORDL-1
        for iEig=1:ORDEig
            for mu=0:min([lam1,lam2,ORDMU-1])

                smallK(lam1+1,lam2+1,mu+1,iEig)=SmallAsympRes(K,iEig,lam1,lam2,mu,3);
                if abs(smallK(lam1+1,lam2+1,mu+1,iEig)) < cutoff
                    out(lam1+1,lam2+1,mu+1,iEig)=smallK(lam1+1,lam2+1,mu+1,iEig);
                else
                    out(lam1+1,lam2+1,mu+1,iEig)=largeK(lam1+1,lam2+1,mu+1,iEig);
                end
            end
        end
    end
end


end

function GLK=SmallAsympRes(K,iEig,lam1,lam2,mu,d)
% calculate the residual using small k asymptot
% iEig: number of eigenvalues, starts from 1
% lam1: spherical harmonic index 1, starts from 0
% lam2: spherical harmonic index 2, starts from 0
% mu: spherical harmonic index 3, starts from 0
l=iEig-1;
Wi=l*(l+d-2);
Res1=1;
Res2=1;

if l>lam1
    for j=lam1:(l-1)
        Wj=j*(j+d-2);
        ajp1=sqrt((j+1-mu)*(j+1+mu+d-3)/(2*j+d)/(2*j+d-2));
        Res1=Res1*ajp1/(Wi-Wj);
    end
    Res1=Res1*(-1i*K)^(l-lam1);
elseif l<lam1
    for j=(l+1):lam1
        Wj=j*(j+d-2);
        ajp1=sqrt((j-mu)*(j+mu+d-3)/(2*j+d-2)/(2*j+d-4));
        Res1=Res1*ajp1/(Wi-Wj);
    end
    Res1=Res1*(-1i*K)^(lam1-l);
end

if l>lam2
    for j=lam2:(l-1)
        Wj=j*(j+d-2);
        ajp1=sqrt((j+1-mu)*(j+1+mu+d-3)/(2*j+d)/(2*j+d-2));
        Res2=Res2*ajp1/(Wi-Wj);
    end
    Res2=Res2*(-1i*K)^(l-lam2);
elseif l<lam2
    for j=(l+1):lam2
        Wj=j*(j+d-2);
        ajp1=sqrt((j-mu)*(j+mu+d-3)/(2*j+d-2)/(2*j+d-4));
        Res2=Res2*ajp1/(Wi-Wj);
    end
    Res2=Res2*(-1i*K)^(lam2-l);
end

% alternatively use following if lam2=0, mu=0;
% for j=0:(l-1)
%     Wj=j*(j+d-2);
%     ajp1=sqrt((j+1)*(j+1+d-3)/(2*j+d)/(2*j+d-2));
%     Res2=Res2*ajp1^2/(Wi-Wj)^2;
% end
% Res2=Res2*K^(2*l)*(-1)^l;
% Res2=sqrt(Res2);

GLK=Res1*Res2;
end


function GLMK=glm(K,EigK,ORDEig,ORDL,NumLayer,ORDMU)
% output GLM(lamda1+1,lamda2+1,mu+1,l+1)
% this function does not contain the small K limit
% K is the magnitude
% EigK: eigenvalues associated with K.  This is a ORDEig x ORDL matrix.
newVersion=1;
if newVersion
MIN=1e-10;

if abs(K)<MIN
    GLMK=zeros(ORDL,ORDL,ORDMU,ORDEig);
    for M=0:(ORDMU-1)
        for L1=M:(ORDL-1)
            GLMK(L1+1,L1+1,M+1,:)=1;
        end
    end
    return
end
if NumLayer-3 < ORDL
    error('make NumLayer bigger')
end
GLMK=zeros(ORDL,ORDL,ORDMU,ORDEig)*NaN; % lam1+1, lam2+1, mu+1, L+1 
for M=0:(ORDMU-1) % M is mu
    % pre calculate a, rows refer to lambda+1
    lamVec=0:NumLayer+1;
    a=sqrt((lamVec-M).*(lamVec+M)./(4*lamVec.^2-1));
    
    for L=M:(ORDEig-1) % may be able to vectorize this loop for speed
        % pre calculate P, rows refer to lambda+1
        P=EigK(L+1,M+1)+lamVec.*(lamVec+1);
        
        % pre calculate \tilde j_\lamda^{\mu (+)}
        % row will refer to \lamda+1
        jp=zeros(NumLayer+1,1)*NaN;
        jp(NumLayer+1)=P(NumLayer+1)/K; % initialize recursion
        for lamPrime=(NumLayer-1):-1:M
            jp(lamPrime+1)=(P(lamPrime+1)/K) + ((a(lamPrime+1+1))^2 / jp(lamPrime+1+1));
        end
        
        % pre calculate \tilde j_\lamda^{\mu (-)} 
        % row will refer to \lamda+1
        jm=zeros(ORDL+2,1)*NaN; % only need jm up to lam0+1
        jm(M+1)=P(M+1)/K; % initialize recursion
        for lamPrime=(M+1):(ORDL+1)
            jm(lamPrime+1)=P(lamPrime+1)/K + (a(lamPrime+1))^2 / jm(lamPrime+1-1);
        end
        
        % pre calculate \partial_p \tilde j_\lamda^{\mu (+)}
        % row will refer to \lamda+1
        djp=zeros(NumLayer+1,1)*NaN;
        djp(NumLayer+1)=1/K;
        for lamPrime=(NumLayer-1):-1:M
            djp(lamPrime+1)=1/K - (a(lamPrime+1+1))^2 * djp(lamPrime+1+1) ...
                            /(jp(lamPrime+1+1)^2);
        end
        
        % pre calculate \tilde j_\lamda^{\mu (-)} 
        % row will refer to \lamda+1
        djm=zeros(ORDL+2,1)*NaN; % only need jm up to lam0+1
        djm(M+1)=1/K; % initialize recursion
        for lamPrime=(M+1):(ORDL+1)
            djm(lamPrime+1)=1/K - (a(lamPrime+1))^2 * djm(lamPrime+1-1) ...
                           / (jm(lamPrime+1-1)^2);
        end
        
        for lam0=M:(ORDL-1)
            if lam0==M % then the first term is zero
                X1prime=( 0 + 1  ...
                          -a(lam0+2)^2*K*djp(lam0+1+1)/((jp(lam0+1+1))^2) );
                if a(lam0+1) ~= 0
                    error('should be zero')
                end
            else
                X1prime=( -a(lam0+1)^2*K*djm(lam0+1-1)/((jm(lam0+1-1))^2) +1  ...
                          -a(lam0+2)^2*K*djp(lam0+1+1)/((jp(lam0+1+1))^2) );
            end

            for lam=M:(ORDL-1)
                if lam==lam0
                    G=1/X1prime;
                elseif lam0<lam
                    % ranges from lam0+1 up to lam
                    X3=prod(jp((lam0+1+1):(lam+1)));
                    G=1i^(lam-lam0)*prod(a((lam0+1+1):(lam+1)))/(X1prime*X3);
                else  % lam0>lam
                    % ranges from lam0 -1 down to and including lam
                    X2=prod(jm((lam+1):(lam0+1-1)));
                    G=1i^(lam0-lam)*prod(a((lam+1+1):(lam0+1-1+1)))/(X1prime*X2);
                    
%                     if isnan(X2)
%                         sprintf('lam0=%d, lam=%d, M=%d, L=%d, K=%g',lam0,lam,M,L,K)
%                         sprintf('X2=%g+%gi, G=%g+%gi, X1prime=%g+%gi',...
%                                 real(X2),imag(X2),real(G),imag(G),real(G),imag(G))
%                         format long
%                         disp('Eigenvalue')
%                         disp(EigK(L+1,M+1))
%                         disp('P')
%                         disp(P(1:10))
%                         disp('jm')
%                         disp(jm)
%                         disp('djm')
%                         disp(djm)
%                         disp('jp')
%                         disp(jp(1:10))
%                         disp('djp')
%                         disp(djp(1:10))
%                         
%                         error('here')
%                     end
                    
                end
                GLMK(lam0+1,lam+1,M+1,L+1)=G;
            end
        end
    end
end


else % end of new version, begging of old

GLMK=zeros(ORDL,ORDL,ORDMU,ORDEig);
MIN=1e-10;

if abs(K)<MIN
    for M=0:(ORDMU-1)
        for L1=M:(ORDL-1)
            GLMK(L1+1,L1+1,M+1,:)=1;
        end
    end
else
 
for M=0:(ORDMU-1)
    AL=zeros(NumLayer,1);
    for iEig=1:ORDEig
        
        WP=zeros(NumLayer,1);
        WM=zeros(NumLayer,1);
        dJp=zeros(NumLayer,1);
        dJm=zeros(NumLayer,1);
        WPPROD=zeros(NumLayer,1);
        
        LP=NumLayer-1;
        AL(NumLayer)=sqrt((LP-M)*(LP+M))/sqrt(4*LP^2-1);
        WP(NumLayer)=1/(EigK(iEig,M+1)+LP*(LP+1)); % Q.J.M. changed 8/3/15
        
        LM=M;
        if EigK(iEig,M+1)+LM*(LM+1)==0 % Q.J.M. changed 8/3/15
            WM(M+1)=0;
        else
            WM(M+1)=1/(EigK(iEig,M+1)+LM*(LM+1)); % Q.J.M. changed 8/3/15
        end
        AL(M+1)=sqrt((LM-M)*(LM+M))/sqrt(4*LM^2-1);
        
        dJp(NumLayer)=1;
        dJm(M+1)=1;
        WPPROD(NumLayer)=1i*K*AL(NumLayer)*WP(NumLayer);
        
        for n=(NumLayer-2):-1:M
            IP=n+1;
            LP=IP-1;
            IM=NumLayer-n+M;
            LM=IM-1;
            
            PL=EigK(iEig,M+1)+LP*(LP+1); % Q.J.M. changed 8/3/15
            AL(IP)=sqrt((LP-M)*(LP+M))/sqrt(4*LP^2-1);
            WP(IP)=1/(PL+WP(IP+1)*(AL(IP+1)*K)^2);
            dJp(IP)=1-dJp(IP+1)*(AL(IP+1)*K*WP(IP+1))^2;
            
            WPPROD(IP)=1i*K*AL(IP)*WP(IP);
            
            PL=EigK(iEig,M+1)+LM*(LM+1); % Q.J.M. changed 8/3/15
            AL(IM)=sqrt((LM-M)*(LM+M))/sqrt(4*LM^2-1);
            if WM(IM-1)==0
                WM(IM)=0;
            else
                WM(IM)=1/(PL+WM(IM-1)*(AL(IM)*K)^2);
            end
            dJm(IM)=1-dJm(IM-1)*(AL(IM)*K*WM(IM-1))^2;
            
        end
        
        for L1=M:(ORDL-1)
            for L2=L1:(ORDL-1)
                L=sort([L1 L2]);
                IM=L(1)+1;
                IP=L(2)+1;
                if L1==M
                    dWL0M=1/dJp(IM);
                else
                    dWL0M=1/(-(AL(IM)*K*WM(IM-1))^2*dJm(IM-1)+dJp(IM));
                end

                if IM==IP
                    GLMK(L1+1,L2+1,M+1,iEig)=dWL0M;
                else
                    GLMK(L1+1,L2+1,M+1,iEig)=dWL0M*prod(WPPROD((IM+1):IP));
                end
                GLMK(L2+1,L1+1,M+1,iEig)=GLMK(L1+1,L2+1,M+1,iEig);               
            end
        end
        
    end
end
end
end % end of old version
end 
