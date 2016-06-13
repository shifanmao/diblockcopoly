function valeq=S3_case2_intOld(FA,R1,R2)
    FB=1-FA;
    MIN=(10^-5)/FA;
    if abs(R1-R2)<MIN
        if abs(R1)<MIN
            valeq=(-1/2).*((-1)+FA).*FA.^2;
        else
            valeq=exp(1).^((-1).*FA.*R1).*(exp(1).^R1+(-1).*exp(1).^(FA.*R1)).*R1.^( ...
              -3).*(1+exp(1).^(FA.*R1).*((-1)+FA.*R1));
        end
    else
        if abs(R1)<MIN
            valeq=exp(1).^((-1).*FA.*R2).*(exp(1).^R2+(-1).*exp(1).^(FA.*R2)).*R2.^( ...
                -3).*((-1)+exp(1).^(FA.*R2)+(-1).*FA.*R2);

        elseif abs(R2)<MIN
            valeq=((-1)+FA).*R1.^(-2).*(1+(-1).*exp(1).^(FA.*R1)+FA.*R1);

        else
            valeq=exp(1).^((-1).*FA.*R2).*(exp(1).^R2+(-1).*exp(1).^(FA.*R2)).*R1.^( ...
              -1).*R2.^(-2).*((-1).*R1+R2).^(-1).*(((-1)+exp(1).^(FA.*R2)).*R1+ ...
              R2+(-1).*exp(1).^(FA.*R1).*R2);

        end
    end
end