function S4=s4wlc(NM,FA,Q1,Q2,Q3,Q4)
%% Function :: s3wlc. Calcultes fourier transform of three-point 
% correlation function of worm-like chain
% Usage: s4wlc(NM,FA,Q1,Q2,Q3,Q4)
% For example:
%   S4(1,1,1,1)=SAAAA
%   S4(1,1,2,1)=SAABA

% parameters for worm-like chain calculations
ORDEig=4;  % maximum number of eigenvalues
ORDL=4;
NumLayer=500;  % number of residual layers
FB=1-FA;
MIN=1e-5;

S4=zeros(2,2,2,2);
% Begin calculation of s4
if sum(power(Q1+Q2+Q3+Q4,2)) <= MIN
    % Evaluate the quantities for s4 calculation
    Q12=Q1+Q2;
    Q13=Q1+Q3;
    Q14=Q1+Q4;
    Q23=Q2+Q3;
    Q24=Q2+Q4;
    Q34=Q3+Q4;

    Q1MAG=sqrt(sum(power(Q1,2)));
    Q2MAG=sqrt(sum(power(Q2,2)));
    Q3MAG=sqrt(sum(power(Q3,2)));
    Q4MAG=sqrt(sum(power(Q4,2)));
    
    R1=Eigenvalues(norm(Q1MAG),ORDEig,ORDL);
    R2=Eigenvalues(norm(Q2MAG),ORDEig,ORDL);
    R3=Eigenvalues(norm(Q3MAG),ORDEig,ORDL);
    R4=Eigenvalues(norm(Q4MAG),ORDEig,ORDL);
    
    Q12MAG=sqrt(sum(power(Q12,2)));
    Q13MAG=sqrt(sum(power(Q13,2)));
    Q14MAG=sqrt(sum(power(Q14,2)));
    Q23MAG=sqrt(sum(power(Q23,2)));
    Q24MAG=sqrt(sum(power(Q24,2)));
    Q34MAG=sqrt(sum(power(Q34,2)));

    R12=Eigenvalues(Q12MAG,ORDEig,ORDL);
    R13=Eigenvalues(Q13MAG,ORDEig,ORDL);
    R14=Eigenvalues(Q14MAG,ORDEig,ORDL);
    R23=Eigenvalues(Q23MAG,ORDEig,ORDL);
    R24=Eigenvalues(Q24MAG,ORDEig,ORDL);
    R34=Eigenvalues(Q34MAG,ORDEig,ORDL);
    
    GL1=Residues(Q1MAG,R1,ORDEig,ORDL,NumLayer,1);
    GL2=Residues(Q2MAG,R2,ORDEig,ORDL,NumLayer,1);
    GL3=Residues(Q3MAG,R3,ORDEig,ORDL,NumLayer,1);
    GL4=Residues(Q4MAG,R4,ORDEig,ORDL,NumLayer,1);
    
    GLM12=Residues(Q12MAG,R12,ORDEig,ORDL,NumLayer,ORDL);
    GLM13=Residues(Q13MAG,R13,ORDEig,ORDL,NumLayer,ORDL);
    GLM14=Residues(Q14MAG,R14,ORDEig,ORDL,NumLayer,ORDL);
    GLM23=Residues(Q23MAG,R23,ORDEig,ORDL,NumLayer,ORDL);
    GLM24=Residues(Q24MAG,R24,ORDEig,ORDL,NumLayer,ORDL);
    GLM34=Residues(Q34MAG,R34,ORDEig,ORDL,NumLayer,ORDL);

    [YLM112,YLM443]=WignerD(Q1,Q12,Q4,ORDL);
    [~,YLM334]=WignerD(Q1,Q12,Q3,ORDL);
    [YLM113,YLM442]=WignerD(Q1,Q13,Q4,ORDL);
    [~,YLM224]=WignerD(Q1,Q13,Q2,ORDL);
    [YLM114,YLM332]=WignerD(Q1,Q14,Q3,ORDL);
    [~,YLM223]=WignerD(Q1,Q14,Q2,ORDL);
    [~,YLM441]=WignerD(Q2,Q23,Q4,ORDL);
    [~,YLM331]=WignerD(Q2,Q24,Q3,ORDL);
    [~,YLM221]=WignerD(Q3,Q34,Q2,ORDL);
    
    R1=NM*R1;
    R2=NM*R2;
    R3=NM*R3;
    R4=NM*R4;
    R12=NM*R12;
    R13=NM*R13;
    R14=NM*R14;
    R23=NM*R23;
    R24=NM*R24;
    R34=NM*R34;
    for M=0:(ORDL-1)
        if mod(ORDEig-M,2)==0
            Lmax=ORDEig-1;
        else
            Lmax=ORDEig-2;
        end
        
        for N1=0:Lmax
            for N2=M:Lmax
                for N3=0:Lmax

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
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R12,R4,FA,FB,GL1,GLM12,GL4,YLM112,YLM443);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R12,R3,FA,FB,GL1,GLM12,GL3,YLM112,YLM334);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R13,R4,FA,FB,GL1,GLM13,GL4,YLM113,YLM442);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R13,R2,FA,FB,GL1,GLM13,GL2,YLM113,YLM224);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R14,R3,FA,FB,GL1,GLM14,GL3,YLM114,YLM332);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R14,R2,FA,FB,GL1,GLM14,GL2,YLM114,YLM223);

                    S4=S4_case1(S4,NM,N1,N2,N3,M,R2,R23,R4,FA,FB,GL2,GLM23,GL4,YLM223,YLM441);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R2,R23,R1,FA,FB,GL2,GLM23,GL1,YLM223,YLM114);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R2,R24,R3,FA,FB,GL2,GLM24,GL3,YLM224,YLM331);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R2,R24,R1,FA,FB,GL2,GLM24,GL1,YLM224,YLM113);

                    S4=S4_case1(S4,NM,N1,N2,N3,M,R3,R34,R2,FA,FB,GL3,GLM34,GL2,YLM334,YLM221);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R3,R34,R1,FA,FB,GL3,GLM34,GL1,YLM334,YLM112);

                    % Case 2: A1=A2=A3~=A4 (SAAAB,SBBBA)
                    S4=S4_case2(S4,NM,N1,N2,N3,M,R1,R12,R4,FA,FB,GL1,GLM12,GL4,YLM112,YLM443);
                    S4=S4_case2(S4,NM,N1,N2,N3,M,R1,R13,R4,FA,FB,GL1,GLM13,GL4,YLM113,YLM442);
                    S4=S4_case2(S4,NM,N1,N2,N3,M,R2,R23,R4,FA,FB,GL2,GLM23,GL4,YLM223,YLM441);

                    % Case 3: A1=A2~=A3=A4 (SAABB)
                    S4=S4_case3(S4,NM,N1,N2,N3,M,R1,R12,R4,FA,GL1,GLM12,GL4,YLM112,YLM443);
                    S4=S4_case3(S4,NM,N1,N2,N3,M,R1,R12,R3,FA,GL1,GLM12,GL3,YLM112,YLM334);
                    S4=S4_case3(S4,NM,N1,N2,N3,M,R2,R12,R4,FA,GL2,GLM12,GL4,YLM221,YLM443);
                    S4=S4_case3(S4,NM,N1,N2,N3,M,R2,R12,R3,FA,GL2,GLM12,GL3,YLM221,YLM334);

                    % Case 4: A1~=A2=A3~=A4 (SABBA)
                    S4=S4_case4(S4,NM,N1,N2,N3,M,R1,R14,R2,FA,GL1,GLM14,GL2,YLM114,YLM223);
                    S4=S4_case4(S4,NM,N1,N2,N3,M,R1,R14,R3,FA,GL1,GLM14,GL3,YLM114,YLM332);
                    S4=S4_case4(S4,NM,N1,N2,N3,M,R4,R14,R2,FA,GL4,GLM14,GL2,YLM441,YLM223);
                    S4=S4_case4(S4,NM,N1,N2,N3,M,R4,R14,R3,FA,GL4,GLM14,GL3,YLM441,YLM332);

                    % Case 5: A1~=A2~=A3~=A4 (SABAB)
                    S4=S4_case5(S4,NM,N1,N2,N3,M,R1,R13,R4,FA,GL1,GLM13,GL4,YLM113,YLM442);
                    S4=S4_case5(S4,NM,N1,N2,N3,M,R1,R13,R2,FA,GL1,GLM13,GL2,YLM113,YLM224);
                    S4=S4_case5(S4,NM,N1,N2,N3,M,R3,R13,R4,FA,GL3,GLM13,GL4,YLM331,YLM442);
                    S4=S4_case5(S4,NM,N1,N2,N3,M,R3,R13,R2,FA,GL3,GLM13,GL2,YLM331,YLM224);
    
                end
            end
        end
    end
    
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
function S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R12,R3,FA,FB,GL1,GLM12,GL3,YLM112,YLM443)
% Case 1: A1=A2=A3=A4
[~,ORDL]=size(R12);
if M==0
    for lam2=M:(ORDL-1)
        for lam3=M:(ORDL-1)
            S4(1,1,1,1)=S4(1,1,1,1)+2*S4_case1_int(FA,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1);
            S4(2,2,2,2)=S4(2,2,2,2)+2*S4_case1_int(FB,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1);        
        end
    end
else
    for lam2=M:(ORDL-1)
        for lam3=M:(ORDL-1)
            S4(1,1,1,1)=S4(1,1,1,1)+4*S4_case1_int(FA,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*real(...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1));
            S4(2,2,2,2)=S4(2,2,2,2)+4*S4_case1_int(FB,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*real(...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1));
        end
    end
end
end
function S4=S4_case2(S4,NM,N1,N2,N3,M,R1,R12,R3,FA,FB,GL1,GLM12,GL3,YLM112,YLM443)
% Case 2: A1=A2=A3~=A4
[~,ORDL]=size(R12);
if M==0
    for lam2=M:(ORDL-1)
        for lam3=M:(ORDL-1)
            S4(1,1,1,2)=S4(1,1,1,2)+2*S4_case2_int(FA,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1);
            S4(2,2,2,1)=S4(2,2,2,1)+2*S4_case2_int(FB,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1);
        end
    end
else
    for lam2=M:(ORDL-1)
        for lam3=M:(ORDL-1)
            S4(1,1,1,2)=S4(1,1,1,2)+4*S4_case2_int(FA,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*real(...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1));
            S4(2,2,2,1)=S4(2,2,2,1)+4*S4_case2_int(FB,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*real(...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1));
        end
    end
end
end
function S4=S4_case3(S4,NM,N1,N2,N3,M,R1,R12,R3,FA,GL1,GLM12,GL3,YLM112,YLM443)
% Case 3: A1=A2~=A3=A4 (SAABB)
[~,ORDL]=size(R12);
if M==0
    for lam2=M:(ORDL-1)
        for lam3=M:(ORDL-1)
            S4(1,1,2,2)=S4(1,1,2,2)+S4_case3_int(FA,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1);
        end
    end
else
    for lam2=M:(ORDL-1)
        for lam3=M:(ORDL-1)
            S4(1,1,2,2)=S4(1,1,2,2)+2*S4_case3_int(FA,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*real(...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1));
        end
    end
end
end
function S4=S4_case4(S4,NM,N1,N2,N3,M,R1,R12,R3,FA,GL1,GLM12,GL3,YLM112,YLM443)
% Case 4: A1~=A2=A3~=A4 (SABBA)
[~,ORDL]=size(R12);
if M==0
    for lam2=M:(ORDL-1)
        for lam3=M:(ORDL-1)
            S4(1,2,2,1)=S4(1,2,2,1)+S4_case3_int(FA,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1);
        end
    end
else
    for lam2=M:(ORDL-1)
        for lam3=M:(ORDL-1)
            S4(1,2,2,1)=S4(1,2,2,1)+2*S4_case3_int(FA,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*real(...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1));
        end
    end
end
end
function S4=S4_case5(S4,NM,N1,N2,N3,M,R1,R12,R3,FA,GL1,GLM12,GL3,YLM112,YLM443)
% Case 5: A1~=A2~=A3~=A4 (SABAB)
[~,ORDL]=size(R12);
if M==0
    for lam2=M:(ORDL-1)
        for lam3=M:(ORDL-1)
            S4(1,2,1,2)=S4(1,2,1,2)+S4_case3_int(FA,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1);
        end
    end
else
    for lam2=M:(ORDL-1)
        for lam3=M:(ORDL-1)
            S4(1,2,1,2)=S4(1,2,1,2)+2*S4_case3_int(FA,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*real(...
                    YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                    *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1));
        end
    end
end
end