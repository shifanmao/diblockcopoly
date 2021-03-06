function valeq=S4_case2_int(FA,R1,R12,R3)
    offset=10^-12;
    if abs(R1)<offset
        R1=offset;
    end
    if abs(R12)<offset
        R12=offset;
    end
    if abs(R3)<offset
        R3=offset;
    end
    
    % on same monomer
    MIN=(10^-1)/FA;
    if max(abs([R1,R12,R3]))<MIN  % 0~R1~R12~R3
        valeq=(-1/6).*((-1)+FA).*FA.^3;
    elseif max(abs([R1,R3]))<MIN   % 0~R1~R3
        valeq=(1/2).*((-1)+FA).*R12.^(-3).*(2+(-2).*exp(1).^(FA.*R12)+FA.*R12.*( ...
                    2+FA.*R12));
    elseif max(abs([R1,R12]))<MIN  % 0~R1~R12
        valeq=(1/2).*exp(1).^((-1).*FA.*R3).*(exp(1).^R3+(-1).*exp(1).^(FA.*R3)) ...
                    .*R3.^(-4).*((-2)+2.*exp(1).^(FA.*R3)+(-1).*FA.*R3.*(2+FA.*R3));
    elseif max(abs([R12,R3]))<MIN % 0~R12~R3
        valeq=(1/2).*((-1)+FA).*R1.^(-3).*(2+(-2).*exp(1).^(FA.*R1)+FA.*R1.*(2+ ...
                    FA.*R1));
    elseif (abs(R1-R12)<MIN && abs(R12-R3)<MIN)
       valeq=(1/2).*exp(1).^((-1).*FA.*R1).*(exp(1).^R1+(-1).*exp(1).^(FA.*R1)) ...
                  .*R1.^(-4).*((-2)+exp(1).^(FA.*R1).*(2+FA.*R1.*((-2)+FA.*R1)));
    elseif abs(R1-R12)<MIN
       valeq=exp(1).^((-1).*FA.*R3).*(exp(1).^R3+(-1).*exp(1).^(FA.*R3)).*R1.^( ...
                  -2).*(R1+(-1).*R3).^(-2).*R3.^(-2).*((-1).*R1.^2+exp(1).^(FA.*R3) ...
                  .*R1.^2+2.*R1.*R3+(-1).*R3.^2+exp(1).^(FA.*R1).*R3.*(R1.*((-2)+ ...
                  FA.*(R1+(-1).*R3))+R3));
    elseif abs(R1-R3)<MIN
       valeq=exp(1).^((-1).*FA.*R1).*(exp(1).^R1+(-1).*exp(1).^(FA.*R1)).*R1.^( ...
                  -3).*(R1+(-1).*R12).^(-2).*R12.^(-1).*(((-1)+exp(1).^(FA.*R1)).* ...
                  R12.^2+R1.^2.*((-1)+exp(1).^(FA.*R12)+exp(1).^(FA.*R1).*FA.*R12)+ ...
                  R1.*R12.*(2+(-1).*exp(1).^(FA.*R1).*(2+FA.*R12)));
    elseif abs(R12-R3)<MIN
       valeq=exp(1).^((-1).*FA.*R3).*(exp(1).^R3+(-1).*exp(1).^(FA.*R3)).*R1.^( ...
                  -1).*(R1+(-1).*R3).^(-2).*R3.^(-3).*((-1).*R1.^2+2.*R1.*R3+(-1).* ...
                  R3.^2+exp(1).^(FA.*R1).*R3.^2+exp(1).^(FA.*R3).*R1.*(R1+(-1).*FA.* ...
                  R1.*R3+R3.*((-2)+FA.*R3)));
    else
       valeq=exp(1).^((-1).*FA.*R3).*(exp(1).^R3+(-1).*exp(1).^(FA.*R3)).*R1.^( ...
                  -1).*(R1+(-1).*R12).^(-1).*R12.^(-1).*(R1+(-1).*R3).^(-1).*(R12+( ...
                  -1).*R3).^(-1).*R3.^(-2).*(((-1)+exp(1).^(FA.*R1)).*R12.*(R12+(-1) ...
                  .*R3).*R3+R1.^2.*(((-1)+exp(1).^(FA.*R3)).*R12+R3+(-1).*exp(1).^( ...
                  FA.*R12).*R3)+R1.*((1+(-1).*exp(1).^(FA.*R3)).*R12.^2+((-1)+exp(1) ...
                  .^(FA.*R12)).*R3.^2));
    end 
end