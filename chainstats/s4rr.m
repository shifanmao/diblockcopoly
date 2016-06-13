function s4=s4rr(N,FA,Q1,Q2,Q3,Q4)

s4=zeros(2,2,2,2);
MIN=1e-5;

if sum(power(Q1+Q2+Q3+Q4,2)) > MIN
    disp('ERROR :: Qs violate translational invariance')
else
    
    for I1=1:N
        for I2=1:N
            for I3=1:N
                for I4=1:N
                  Qtot=norm(I1*Q1+I2*Q2+I3*Q3+I4*Q4);
                  pgc=zeros(2,2,2,2);

                  for A1=1:2
                      for A2=1:2
                          for A3=1:2
                            for A4=1:2
                              if Qtot<1e-5
                                  pgc(A1,A2,A3,A4)=1;
                              else
                                  pgc(A1,A2,A3,A4)=sin(Qtot)/(Qtot);
                              end
                              s4(A1,A2,A3,A4)=s4(A1,A2,A3,A4)+...
                                  pa1a2a3a4(A1,A2,A3,A4,FA,N,I1,I2,I3,I4)*pgc(A1,A2,A3,A4);
                            end
                          end
                      end
                  end

                end
            end
        end
    end

end

end
 
function [PA1A2A3A4]=pa1a2a3a4(A1,A2,A3,A4,FA,N,I1,I2,I3,I4)
    IND1=(I1-1)/N;
    IND2=(I2-1)/N;
    IND3=(I3-1)/N;
    IND4=(I4-1)/N;
    PA1A2A3A4=0;
    
    if ( ((IND1<FA)&&(A1==1) || (IND1>=FA)&&(A1==2)) &&...
         ((IND2<FA)&&(A2==1) || (IND2>=FA)&&(A2==2)) &&...
         ((IND3<FA)&&(A3==1) || (IND3>=FA)&&(A3==2)) &&...
         ((IND4<FA)&&(A4==1) || (IND4>=FA)&&(A4==2)))
        PA1A2A3A4=1;
    end
end
