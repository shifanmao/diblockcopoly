function valeq=S3_case1_int(FA,R1,R2)
    MIN=(10^-2)/FA;
    if abs(R1-R2)<MIN
        if abs(R1)<MIN
            valeq=(1/6).*FA.^3;
        else
            valeq=R1.^(-3).*(2+FA.*R1+exp(1).^(FA.*R1).*((-2)+FA.*R1));
        end
    else
        if abs(R1)<MIN
            valeq=(-1/2).*R2.^(-3).*(2+(-2).*exp(1).^(FA.*R2)+2.*FA.*R2+FA.^2.*R2.^2);
        elseif abs(R2)<MIN
            valeq=(-1/2).*R1.^(-3).*(2+(-2).*exp(1).^(FA.*R1)+2.*FA.*R1+FA.^2.*R1.^2);
        else
            valeq=R1.^(-2).*(R1+(-1).*R2).^(-1).*R2.^(-2).*(((-1)+exp(1).^(FA.*R1)) ...
              .*R2.^2+(-1).*FA.*R1.*R2.^2+R1.^2.*(1+(-1).*exp(1).^(FA.*R2)+FA.*R2));
        end
    end
end