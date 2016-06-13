function S3=s3wlc(NM,FA,Q1,Q2,Q3)
%% Function :: s3wlc. Calcultes fourier transform of three-point 
% correlation function of worm-like chain
% Usage: S3=s3wlc(NM,FA,Q1,Q2,Q3,ORDEig,ORDL,ResLayer)
% For example:
%   S3(1,1,1)=SAAA
%   S3(1,1,2)=SAAB

% parameters for worm-like chain calculations
ORDEig=20;  % maximum number of eigenvalues
ORDL=20;
ResLayer=500;  % number of residual layers
FB=1-FA;

NR=ORDEig;
MIN=1e-5;

S3=zeros(2,2,2);
% Begin calculation of s3
if sum(power(Q1+Q2+Q3,2)) > MIN
    
    disp(['sum(Q)=',num2str(sum(power(Q1+Q2+Q3,2)))])
    error('Wavevectors must add up to zero from translational invariance')
    
elseif sum(power(Q3,2)) < MIN
    
    % reduces to two-point correlation function
    Q1MAG=sqrt(sum(power(Q1,2)));
    S2 = s2wlc(NM,FA,Q1MAG);
    S3(1,1,1) = S2(1,1)*FA*NM;
    S3(1,1,2) = S2(1,1)*FB*NM;
    S3(1,2,1) = S2(1,2)*FA*NM;
    S3(1,2,2) = S2(1,2)*FB*NM;
    S3(2,1,1) = S2(2,1)*FA*NM;
    S3(2,1,2) = S2(2,1)*FB*NM;
    S3(2,2,1) = S2(2,2)*FA*NM;
    S3(2,2,2) = S2(2,2)*FB*NM;
    
else

    % Evaluate the quantities for s3 calculation
    Q1MAG=sqrt(sum(power(Q1,2)));
    Q2MAG=sqrt(sum(power(Q2,2)));
    Q3MAG=sqrt(sum(power(Q3,2)));

    % calculate the eigenvalues
    R1=Eigenvalues(Q1MAG,NR,1);
    R2=Eigenvalues(Q2MAG,NR,1);
    R3=Eigenvalues(Q3MAG,NR,1);

    % get the residues for all roots of each k(j)
    GL1=Residues(Q1MAG,R1,ORDEig,ORDL,ResLayer,1);
    GL2=Residues(Q2MAG,R2,ORDEig,ORDL,ResLayer,1);
    GL3=Residues(Q3MAG,R3,ORDEig,ORDL,ResLayer,1);

    % but only need one -> m=0
    GL1=GL1(:,1,1,:);
    GL2=GL2(:,1,1,:);
    GL3=GL3(:,1,1,:);

    % unit vectors in direction of wavevectors
    EQ1=Q1/Q1MAG;
    EQ2=Q2/Q2MAG;
    EQ3=Q3/Q3MAG;

    % angles between wavevectors
    RHO12=sum(EQ1.*EQ2);
    RHO23=sum(EQ2.*EQ3);
    RHO13=sum(EQ1.*EQ3);

    % legendre polynomial representations of angles
    PL12=legendrep(-RHO12,ORDL);
    PL23=legendrep(-RHO23,ORDL);
    PL13=legendrep(-RHO13,ORDL);

    R1=NM*R1;
    R2=NM*R2;
    R3=NM*R3;
    for N1=1:NR
        for N2=1:NR
            % Case 1: A1=A2=A3 (SAAA)
            S3=S3_case1(S3,NM,N1,N2,R1,R2,GL1,GL2,PL12,FA,FB);
            S3=S3_case1(S3,NM,N1,N2,R1,R3,GL1,GL3,PL13,FA,FB);
            S3=S3_case1(S3,NM,N1,N2,R2,R3,GL2,GL3,PL23,FA,FB);
            S3=S3_case1(S3,NM,N1,N2,R2,R1,GL2,GL1,PL12,FA,FB);
            S3=S3_case1(S3,NM,N1,N2,R3,R1,GL3,GL1,PL13,FA,FB);
            S3=S3_case1(S3,NM,N1,N2,R3,R2,GL3,GL2,PL23,FA,FB);

            % Case 2: A1=A2~=A3 (SAAB)
            S3=S3_case2(S3,NM,N1,N2,R1,R3,GL1,GL3,PL13,FA,FB);
            S3=S3_case2(S3,NM,N1,N2,R2,R3,GL2,GL3,PL23,FA,FB);

            % Case 3: A1~=A2=A3 (SABB)
            S3=S3_case3(S3,NM,N1,N2,R1,R2,GL1,GL2,PL12,FA,FB);
            S3=S3_case3(S3,NM,N1,N2,R1,R3,GL1,GL3,PL13,FA,FB);
        end
    end

    S3(1,2,1)=S3(2,1,1);
    S3(2,1,2)=S3(1,2,2);
    S3(isnan(S3))=0;

end
end
function S3=S3_case1(S3,NM,N1,N2,R1,R2,GL1,GL2,PL12,FA,FB)
% Case 1: A1=A2=A3
    S3(1,1,1)=S3(1,1,1)+S3_case1_int(FA,R1(N1),R2(N2))*NM^3*...
                                    sum(PL12.*GL1(:,N1).*GL2(:,N2));
    S3(2,2,2)=S3(2,2,2)+S3_case1_int(FB,R1(N1),R2(N2))*NM^3*...
                                    sum(PL12.*GL1(:,N1).*GL2(:,N2));
end
function S3=S3_case2(S3,NM,N1,N2,R1,R2,GL1,GL2,PL12,FA,FB)
% Case 2: A1=A2~=A3
    S3(1,1,2)=S3(1,1,2)+S3_case2_int(FA,R1(N1),R2(N2))*NM^3*...
                                    sum(PL12.*GL1(:,N1).*GL2(:,N2));
    S3(2,2,1)=S3(2,2,1)+S3_case2_int(FB,R1(N1),R2(N2))*NM^3*...
                                    sum(PL12.*GL1(:,N1).*GL2(:,N2));
end
function S3=S3_case3(S3,NM,N1,N2,R1,R2,GL1,GL2,PL12,FA,FB)
% Case 3: A1~=A2=A3
    S3(2,1,1)=S3(2,1,1)+S3_case3_int(FB,R1(N1),R2(N2))*NM^3*...
                                    sum(PL12.*GL1(:,N1).*GL2(:,N2));
    S3(1,2,2)=S3(1,2,2)+S3_case3_int(FA,R1(N1),R2(N2))*NM^3*...
                                    sum(PL12.*GL1(:,N1).*GL2(:,N2));
end