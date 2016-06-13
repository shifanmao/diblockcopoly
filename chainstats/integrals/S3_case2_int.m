function valeq=S3_case2_int(FA,R1,R2)
    FB=1-FA;
    MIN=(10^-5);
    if abs(R1-R2)<MIN
        if abs(R1)<MIN
            % valeq=(-1/2).*((-1)+FA).*FA.^2;
            valeq=FA^2*FB*(3*R1+2*FA*R1-FA*R1+6)/12;
        else
%             valeq=exp(1).^((-1).*FA.*R1).*(exp(1).^R1+(-1).*exp(1).^(FA.*R1)).*R1.^( ...
%               -3).*(1+exp(1).^(FA.*R1).*((-1)+FA.*R1));
            valeq=(expl(3,R1*FB)+expl(3,FA*R1)-expl(3,R1)-FA*R1*expl(2,FA*R1)+FA*R1*expl(2,R1))/(R1^3);
        end
    else
        if abs(R1)<MIN
%             valeq=exp(1).^((-1).*FA.*R2).*(exp(1).^R2+(-1).*exp(1).^(FA.*R2)).*R2.^( ...
%                 -3).*((-1)+exp(1).^(FA.*R2)+(-1).*FA.*R2);
            valeq=-(expl(3,R2*FB)+expl(3,FA*R2)-expl(3,R2)+FA*R2*expl(2,R2*FB))/(R2^3);
        elseif abs(R2)<MIN
%             valeq=((-1)+FA).*R1.^(-2).*(1+(-1).*exp(1).^(FA.*R1)+FA.*R1);
            valeq=FB*expl(2,FA*R1)/R1^2;
        else
%             valeq=exp(1).^((-1).*FA.*R2).*(exp(1).^R2+(-1).*exp(1).^(FA.*R2)).*R1.^( ...
%               -1).*R2.^(-2).*((-1).*R1+R2).^(-1).*(((-1)+exp(1).^(FA.*R2)).*R1+ ...
%               R2+(-1).*exp(1).^(FA.*R1).*R2);
            valeq=(exp(R2*FB)-1)*(R1*expl(1,FA*R2)-R2*expl(1,FA*R1))/(R1*R2^2*(R2-R1));
        end
    end
end