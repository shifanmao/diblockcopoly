function S3=s3gc(N,FA,Q1,Q2,Q3)
%S3 is a three point correlation function
%For example:
%   S3(1,1,1)=SAAA
%   S3(1,1,2)=SAAB

MIN=1e-5;
d=3;    %dimension
FB=1-FA;

% Evaluate the quantities for s3 calculation
Q1MAG=sqrt(sum(power(Q1,2)));
Q2MAG=sqrt(sum(power(Q2,2)));
Q3MAG=sqrt(sum(power(Q3,2)));

S3=zeros(2,2,2);
if sum(power(Q1+Q2+Q3,2)) > MIN
    error('Qs must add up to zero from translationsl invariance')
else
    R1=-N*Q1MAG*Q1MAG/(2*d);
    R2=-N*Q2MAG*Q2MAG/(2*d);
    R3=-N*Q3MAG*Q3MAG/(2*d);

    % Case 1: A1=A2=A3 (SAAA)
    S3=S3_case1(S3,N,R1,R2,FA,FB);
    S3=S3_case1(S3,N,R1,R3,FA,FB);
    S3=S3_case1(S3,N,R2,R3,FA,FB);
    S3=S3_case1(S3,N,R2,R1,FA,FB);
    S3=S3_case1(S3,N,R3,R1,FA,FB);
    S3=S3_case1(S3,N,R3,R2,FA,FB);

    % Case 2: A1=A2~=A3 (SAAB)
    S3=S3_case2(S3,N,R1,R3,FA,FB);
    S3=S3_case2(S3,N,R2,R3,FA,FB);

    % Case 3: A1~=A2=A3 (SABB)
    S3=S3_case3(S3,N,R1,R2,FA,FB);
    S3=S3_case3(S3,N,R1,R3,FA,FB);

    S3(1,2,1)=S3(2,1,1);
    S3(2,1,2)=S3(1,2,2);
    
    S3(isnan(S3))=0;
end
end
function S3=S3_case1(S3,N,R1,R2,FA,FB)
% Case 1: A1=A2=A3
    S3(1,1,1)=S3(1,1,1)+S3_case1_int(FA,R1,R2)*N^3;
    S3(2,2,2)=S3(2,2,2)+S3_case1_int(FB,R1,R2)*N^3;
end
function S3=S3_case2(S3,N,R1,R2,FA,FB)
% Case 2: A1=A2~=A3
    S3(1,1,2)=S3(1,1,2)+S3_case2_int(FA,R1,R2)*N^3;
    S3(2,2,1)=S3(2,2,1)+S3_case2_int(FB,R1,R2)*N^3;
end
function S3=S3_case3(S3,N,R1,R2,FA,FB)
% Case 3: A1~=A2=A3
    S3(2,1,1)=S3(2,1,1)+S3_case3_int(FB,R1,R2)*N^3;
    S3(1,2,2)=S3(1,2,2)+S3_case3_int(FA,R1,R2)*N^3;
end