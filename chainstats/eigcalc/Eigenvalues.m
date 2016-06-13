function Eig=Eigenvalues(K,ORDEig,ORDL)
% Calculates the eigenvalues for the inverse laplace transform
% Returns a ORDEig x ORDL matrix
% Index of output: (l+1,M+1)
% Choose cutoff manually based on previous compairison
% You may not need all output values, i.e. if you know M=0

cutoff = 800; % This appears to be a pretty good cutoff

if K<cutoff
    Eig=IntermediateKEigenValues(K,ORDEig,ORDL);
else
    Eig=LargeKEigenValues(K,ORDEig,ORDL,3);
end

end

function Eig=IntermediateKEigenValues(k,ORDEig,ORDL)
% essentially does the same thing as MatRootsLM
% output size: ORDEig x ORDL
% output indices: (l+1,m+1)

N=2*(ceil(ORDEig/2));
Eig=zeros(N,ORDL);  % output no-longer square

for M=0:(ORDL-1)
    % use matrix method for intermediate and small k regime
    n=4*N;
    E=zeros(n,n);
    for I=1:n
        L=I+M;
        if k<=1
            a=complex(0,-k*sqrt((L-M)*(L+M))/sqrt(4*L^2-1));
            if I>1 
                b=complex(0,-k*sqrt((L-1-M)*(L-1+M))/sqrt(4*(L-1)^2-1)); 
            end
            if I==1
                E(I,1:2)=[L*(L-1),a];
            elseif I==n
                E(I,n-1:n)=[b,L*(L-1)];
            else
                E(I,I-1:I+1)=[b,L*(L-1),a];
            end
        else
            a=complex(0,-sqrt((L-M)*(L+M))/sqrt(4*L^2-1));
            if I>1 
                b=complex(0,-sqrt((L-1-M)*(L-1+M))/sqrt(4*(L-1)^2-1)); 
            end
            if I==1
                E(I,1:2)=[L*(L-1)/k,a];
            elseif I==n
                E(I,n-1:n)=[b,L*(L-1)/k];
            else
                E(I,I-1:I+1)=[b,L*(L-1)/k,a];
            end
        end
    end  
    TempMat=eig(E);
    [junk,index]=sort(real(TempMat));
    TempMat=TempMat(index);
    if N-M<1
        Eig(:,M+1)=zeros(N,1)*NaN;
    elseif k<=1
        Eig(:,M+1)=-TempMat(1:N); 
    else
        Eig(:,M+1)=-TempMat(1:N)*k; 
    end
end

% Inforce choice of how to order l values, choice is somewhat arbitrary
% This used to be done outside of MatRoots
for M=0:ORDL-1
    for I=1:2:ORDEig
        Eig(I,M+1)=real(Eig(I,M+1))+1i*abs(imag(Eig(I,M+1)));
        Eig(I+1,M+1)=real(Eig(I+1,M+1))-1i*abs(imag(Eig(I+1,M+1)));    
    end
end

% delete extra row if it exists
if N>ORDEig
    Eig(end,:)=[];
end

% shift to higher L
for M=0:(ORDL-1)
    if ORDEig-M<1
        Eig(:,M+1)=zeros(ORDEig,1)*NaN;
    else  
        Eig(:,M+1)=[zeros(M,1)*NaN;Eig(1:(ORDEig-M),M+1)]; % Q.J.M. added this 8/16/15
    end
end


end

function Eig=LargeKEigenValues(k,ORDEig,ORDL,d)
% Returs an ORDEig x ORDL  where entries are for l+1,mu+1
% where l is the order of the eigenvalue
% and mu is is the azimuthal quantum number (i.e. m)

Eig=zeros(ORDEig,ORDL)*NaN;
mu=0:(ORDL-1);

for I=1:floor(ORDEig/2)
    l=I-1;
    alpha=1/sqrt(8*k);
    Eig(2*l+1,:)=1i*k-mu.*(mu+d-2)-Epsilon(l,d,alpha,mu); % Q.J.M. changes this line 8/1/15
    Eig(2*l+2,:)=conj(Eig(2*l+1,:));
end

% Q.J.M. added this loop 8/16/15
for M=0:(ORDL-1)
    if ORDEig-M<1
        Eig(:,M+1)=zeros(ORDEig,1)*NaN;
    else  
        Eig(:,M+1)=[zeros(M,1)*NaN;Eig(1:(ORDEig-M),M+1)];
    end
end
    
end

function value=Epsilon(l,d,alpha,mu)
% For use in calculating Large K EigenValues
% eigenvalues using large k asymptotic expansion
% generates epsilon^{s}_r (see the paper)

I=complex(0,1);
beta=-sqrt(2)/4*(1+I);
m=mu+(d-3)/2;  % Q.J.M. changes this line 8/1/15
n=2*l+m+1;  % also know as s

epsilon_0=(-1/2/beta)^(-1)*(n/2);
epsilon_1=(-1/2/beta)^( 0)*(-1/8*(n.^2+3-3*m.^2)-m.*(m+1));
epsilon_2=(-1/2/beta)^( 1)*(-1/2^5*n.*(n.^2+3-9*m.^2));
epsilon_3=(-1/2/beta)^( 2)*(-1/2^8*(5*n.^4+34*n.^2+9)-(102*n.^2+42).*m.^2+33*m.^4);
epsilon_4=(-1/2/beta)^( 3)*(-1/2^11*n.*(33*n.^4+410*n.^2+405)-(1230*n.^2+1722).*m.^2+813*m.^4);
epsilon_5=(-1/2/beta)^( 4)*(-1/2^12*9*(7*n.^6+140*n.^4+327*n.^2+54-(420*n.^4+1350*n.^2+286).*m.^2+(495*n.^2+314).*m.^4-82*m.^6));

value=epsilon_0/alpha+epsilon_1+epsilon_2*alpha+epsilon_3*alpha^2+...
      epsilon_4*alpha^3+epsilon_5*alpha^4;
end