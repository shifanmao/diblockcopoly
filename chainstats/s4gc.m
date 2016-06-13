function S4=s4gc(N,FA,Q1,Q2,Q3,Q4)
%S4 is a four point correlation function
%For example:
%   S4(1,1,1,1)=SAAAA
%   S4(1,1,2,1)=SAABA

MIN=1e-5;
d=3;    %dimension
FB=1-FA;

% Begin calculation of s4
S4=zeros(2,2,2,2);
if sum(power(Q1+Q2+Q3+Q4,2)) > MIN
    error('Qs must add up to zero from translationsl invariance')
else
    % Evaluate the quantities for s4 calculation
    Q1MAG=sqrt(sum(power(Q1,2)));
    Q2MAG=sqrt(sum(power(Q2,2)));
    Q3MAG=sqrt(sum(power(Q3,2)));
    Q4MAG=sqrt(sum(power(Q4,2)));
    
    R1=-N*Q1MAG*Q1MAG/(2*d);
    R2=-N*Q2MAG*Q2MAG/(2*d);
    R3=-N*Q3MAG*Q3MAG/(2*d);
    R4=-N*Q4MAG*Q4MAG/(2*d);
    
    Q12MAG=sqrt(sum(power(Q1+Q2,2)));
    Q13MAG=sqrt(sum(power(Q1+Q3,2)));
    Q14MAG=sqrt(sum(power(Q1+Q4,2)));
    Q23MAG=sqrt(sum(power(Q2+Q3,2)));
    Q24MAG=sqrt(sum(power(Q2+Q4,2)));
    Q34MAG=sqrt(sum(power(Q3+Q4,2)));

    R12=-N*Q12MAG*Q12MAG/(2*d);
    R13=-N*Q13MAG*Q13MAG/(2*d);
    R14=-N*Q14MAG*Q14MAG/(2*d);
    R23=-N*Q23MAG*Q23MAG/(2*d);
    R24=-N*Q24MAG*Q24MAG/(2*d);
    R34=-N*Q34MAG*Q34MAG/(2*d);
         
    % PERMUTATIONS
    % 1 12 4 112 443
    % 1 12 3 112 334
    % 1 13 4 113 442
    % 1 13 2 113 224
    % 1 14 3 114 332
    % 1 14 2 114 223
    % 2 23 4 223 441
    % 2 23 1 223 114
    % 2 24 3 224 331
    % 2 24 1 224 113
    % 3 34 2 334 221
    % 3 34 1 334 112
    
    % Case 1: A1=A2=A3=A4 (SAAAA,SBBBB)
    S4=S4_case1(S4,N,R1,R12,R3,FA,FB);
    S4=S4_case1(S4,N,R1,R12,R4,FA,FB);
    S4=S4_case1(S4,N,R1,R13,R2,FA,FB);
    S4=S4_case1(S4,N,R1,R13,R4,FA,FB);
    S4=S4_case1(S4,N,R1,R14,R2,FA,FB);
    S4=S4_case1(S4,N,R1,R14,R3,FA,FB);

    S4=S4_case1(S4,N,R2,R23,R1,FA,FB);
    S4=S4_case1(S4,N,R2,R23,R4,FA,FB);
    S4=S4_case1(S4,N,R2,R24,R1,FA,FB);
    S4=S4_case1(S4,N,R2,R24,R3,FA,FB);

    S4=S4_case1(S4,N,R3,R34,R1,FA,FB);
    S4=S4_case1(S4,N,R3,R34,R2,FA,FB);
    
    % Case 2: A1=A2=A3~=A4 (SAAAB,SBBBA)
    S4=S4_case2(S4,N,R1,R12,R4,FA,FB);
    S4=S4_case2(S4,N,R1,R13,R4,FA,FB);
    S4=S4_case2(S4,N,R2,R23,R4,FA,FB);
    
    % Case 3: A1=A2~=A3=A4 (SAABB)
    S4=S4_case3(S4,N,R1,R12,R4,FA);
    S4=S4_case3(S4,N,R1,R12,R3,FA);
    S4=S4_case3(S4,N,R2,R12,R4,FA);
    S4=S4_case3(S4,N,R2,R12,R3,FA);
    
    % Case 4: A1~=A2=A3~=A4 (SABBA)
    S4=S4_case4(S4,N,R1,R14,R2,FA);
    S4=S4_case4(S4,N,R1,R14,R3,FA);
    S4=S4_case4(S4,N,R4,R14,R2,FA);
    S4=S4_case4(S4,N,R4,R14,R3,FA);
    
    % Case 5: A1~=A2~=A3~=A4 (SABAB)
    S4=S4_case5(S4,N,R1,R13,R4,FA);
    S4=S4_case5(S4,N,R1,R13,R2,FA);
    S4=S4_case5(S4,N,R3,R13,R4,FA);
    S4=S4_case5(S4,N,R3,R13,R2,FA);
    
    S4(1,1,2,1)=S4(1,1,1,2);
    S4(1,2,1,1)=S4(1,1,1,2);
    S4(2,1,1,1)=S4(1,1,1,2);
    
    S4(2,2,1,1)=S4(1,1,2,2);
    S4(2,1,1,2)=S4(1,2,2,1);
    S4(2,1,2,1)=S4(1,2,1,2);
    
    S4(2,2,1,2)=S4(2,2,2,1);
    S4(2,1,2,2)=S4(2,2,2,1);
    S4(1,2,2,2)=S4(2,2,2,1);
end
end
function S4=S4_case1(S4,N,R1,R12,R3,FA,FB)
% Case 1: A1=A2=A3=A4
    S4(1,1,1,1)=S4(1,1,1,1)+2*S4_case1_int(R1,R12,R3,FA)*N^4;
    S4(2,2,2,2)=S4(2,2,2,2)+2*S4_case1_int(R1,R12,R3,FB)*N^4;
end
function S4=S4_case2(S4,N,R1,R12,R3,FA,FB)
% Case 2: A1=A2=A3~=A4
    S4(1,1,1,2)=S4(1,1,1,2)+2*S4_case2_int(FA,R1,R12,R3)*N^4;
    S4(2,2,2,1)=S4(2,2,2,1)+2*S4_case2_int(FB,R1,R12,R3)*N^4;
end
function S4=S4_case3(S4,N,R1,R12,R3,FA)
% Case 3: A1=A2~=A3=A4 (SAABB)
    S4(1,1,2,2)=S4(1,1,2,2)+S4_case3_int(FA,R1,R12,R3)*N^4;
end
function S4=S4_case4(S4,N,R1,R12,R3,FA)
% Case 4: A1~=A2=A3~=A4 (SABBA)
    S4(1,2,2,1)=S4(1,2,2,1)+S4_case3_int(FA,R1,R12,R3)*N^4;
end
function S4=S4_case5(S4,N,R1,R12,R3,FA)
% Case 5: A1~=A2~=A3~=A4 (SABAB)
    S4(1,2,1,2)=S4(1,2,1,2)+S4_case3_int(FA,R1,R12,R3)*N^4;
end

function valeq=S4_case1_int(E1,E12,E3,N)
MIN=(10^-3)/N;
vec=[E1,E12,E3];
nzeros=sum(abs(vec)<MIN);
if nzeros==3
    valeq=N^4*(N*(E1+E12+E3)+5)/120;
elseif nzeros==2
    [~,index]=max(abs(vec));
    others=vec; others(index)=[];
    valeq=chicken(others(1),others(2),vec(index),N);
elseif nzeros==1
    [~,index]=min(abs(vec));
    others=vec; others(index)=[];
    if abs(others(1)-others(2))<MIN
        val=0.5*(others(1)+others(2));
        valeq=(-6*expl(3,N*val)+2*N*val*expl(2,N*val))/(2*val^4);
    else
        valeq=(expl(4,N*others(2))/(others(2)^3)...
            -expl(4,N*others(1))/(others(1)^3))...
            /(others(2)-others(1));
    end
else
    dif=abs([E12-E3,E3-E1,E1-E12]);
    if max(dif)<MIN
        val=mean([E1,E12,E3]);
        valeq=(1/2)*val^(-4)*((-2)*(3+val*N)+exp(val*N)*(6+val*N*(( ...
            -4)+val*N)));
    elseif min(dif)>MIN
        valeq=E1.^(-2).*(E1+(-1).*E12).^(-1).*E12.^(-2).*(E1+(-1).*E3).^(-1).*( ...
          E12+(-1).*E3).^(-1).*E3.^(-2).*(((-1)+exp(1).^(E1.*N)).*E12.^2.*( ...
          E12+(-1).*E3).*E3.^2+E1.*E12.^2.*E3.^2.*((-1).*E12+E3).*N+E1.^3.* ...
          ((-1).*((-1)+exp(1).^(E12.*N)).*E3.^2+E12.*E3.^2.*N+E12.^2.*(( ...
          -1)+exp(1).^(E3.*N)+(-1).*E3.*N))+E1.^2.*(((-1)+exp(1).^(E12.* ...
          N)).*E3.^3+(-1).*E12.*E3.^3.*N+E12.^3.*(1+(-1).*exp(1).^(E3.*N) ...
          +E3.*N)));
    else
        [~,index]=min(dif);
        val=vec(index);
        others=vec; others(index)=[];
        oval=0.5*sum(others);
        valeq=((-3*oval*expl(4,N*oval)+N*oval^2*expl(3,N*oval))*val^2+...
          (2*expl(4,N*oval)-N*oval*expl(3,N*oval))*val^3+...
          oval^3*expl(4,N*val))/(oval^3*val^2*(oval-val)^2);       
        
    end
    
end
end