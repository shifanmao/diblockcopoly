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
    S4(1,1,1,1)=S4(1,1,1,1)+2*S4_case1_int(FA,R1,R12,R3)*N^4;
    S4(2,2,2,2)=S4(2,2,2,2)+2*S4_case1_int(FB,R1,R12,R3)*N^4;
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
