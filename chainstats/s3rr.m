function s3=s3rr(N,FA,Q1,Q2,Q3)

s3=zeros(2,2,2);
MIN=1e-5;

if sum(power(Q1+Q2+Q3,2)) > MIN
    disp('ERROR :: Qs violate translational invariance')
else

    for I1=1:N
        for I2=1:N
            for I3=1:N
              Qtot=norm(I1*Q1+I2*Q2+I3*Q3);
              pgc=zeros(2,2,2);

                  for A1=1:2
                      for A2=1:2
                          for A3=1:2
                          if Qtot<1e-5
                              pgc(A1,A2,A3)=1;
                          else
                              pgc(A1,A2,A3)=sin(Qtot)/(Qtot);
                          end
                              s3(A1,A2,A3)=s3(A1,A2,A3)+pa1a2a3(A1,A2,A3,FA,N,I1,I2,I3)*pgc(A1,A2,A3);
                          end
                      end
                  end

            end
        end
    end

end

end
 
function [PA1A2A3]=pa1a2a3(A1,A2,A3,FA,N,I1,I2,I3)
    IND1=(I1-1)/N;
    IND2=(I2-1)/N;
    IND3=(I3-1)/N;
    PA1A2A3=0;
    
    if ( ((IND1<FA)&&(A1==1) || (IND1>=FA)&&(A1==2)) &&...
         ((IND2<FA)&&(A2==1) || (IND2>=FA)&&(A2==2)) &&...
         ((IND3<FA)&&(A3==1) || (IND3>=FA)&&(A3==2)) )
        PA1A2A3=1;
    end
end
