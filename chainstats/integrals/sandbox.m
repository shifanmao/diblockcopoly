clc
FA=0.7;
FB=1-FA;
R1=1.347;
R12=1.1234;
R3=1.4854;
R2=R12;

% ---- S4_case2_int ----
disp(' R1=0, R12=R3')
(2*(expl(1,R3)*expl(3,-FA*R3)+2*coshl(4,FA*R3))...
    +FA*R3*(expl(1,R3)*expl(2,-FA*R3)-2*sinhl(3,FA*R3)))/(R3^4)

disp('R12=0, R1=R3')
(2*(expl(1,R3)*expl(3,-FA*R3)+2*coshl(4,FA*R3))...
    +FA*R3*(expl(1,R3)*expl(2,-FA*R3)-2*sinhl(3,FA*R3)))/(R3^4)

disp('R3=0, R1=R12')
(1-FA)*(-2*expl(3,FA*R1)+FA*R1*expl(2,FA*R1))/(R1^3)

disp('R1=0')
(exp(FB*R3)-1)*(R12^2*expl(3,FA*R3)-R3^2*expl(3,FA*R12))/(R3^3*R12^2*(R3-R12))

disp('R12=0')
(exp(FB*R3)-1)*(R1^2*expl(3,FA*R3)-R3^2*expl(3,FA*R1))/(R3^3*R1^2*(R3-R1))

disp('R3=0')
FB*(R1^2*expl(3,FA*R12)-R12^2*expl(3,FA*R1))/(R1^2*R12^2*(R12-R1))

disp('R1=R12=R3=0')
FA^2*FB*(3*R12+2*FA*R1-FA*R12+6)/12

% ---- S4_case3_int ----

disp('R1=0, R12=R3')
(-FA*FB*expl(2,R3*FB)*R3^2+...
 (FB*expl(3,R3)+(2*FA-1)*expl(3,R3*FB))*R3+...
 expl(4,FA*R3)-expl(4,R3)+expl(4,R3*FB) )/(R3^4)

disp('R3=0, R1=R12')
(-FB*FA*expl(2,R1*FA)*R1^2+...
 (FA*expl(3,R1)+(2*FB-1)*expl(3,R1*FA))*R1+...
 expl(4,FB*R1)-expl(4,R1)+expl(4,R1*FA) )/(R1^4)

disp('R12=0, R1=R3')
( R3*(-FB*expl(3,FA*R3)-FA*expl(3,R3*FB))+...
    expl(4,R3)-expl(4,FA*R3)-expl(4,R3*FB) )/(R3^4)

disp('R1=0')
( R3*expl(2,R12*FB)-R12*expl(2,R3*FB) )*expl(2,FA*R12)/(R3*R12^3*(R12-R3))

disp('R3=0')
( R1*expl(2,R12*FA)-R12*expl(2,R1*FA) )*expl(2,FB*R12)/(R1*R12^3*(R12-R1))

disp('R12=0')
expl(2,FA*R1)*expl(2,R3*FB)/(R1^2*R3^2)

disp('R1=R12=R3=0')
FA^2*FB^2*(R3+E12+FA*E1-FA*E3+3)/12

% ---- S2_case2_int -----

disp('S2_case2_int')
(exp(R12*FB)-1)*(R1*expl(1,FA*R12)-R12*expl(1,FA*R1))/(R1*R12^2*(R12-R1))

disp('R1=R2')
(expl(3,R1*FB)+expl(3,FA*R1)-expl(3,R1)-FA*R1*expl(2,FA*R1)+FA*R1*expl(2,R1))/(R1^3)

disp('R1=0')
-(expl(3,R2*FB)+expl(3,FA*R2)-expl(3,R2)+FA*R2*expl(2,R2*FB))/(R2^3)

disp('R2=0')
FB*expl(2,FA*R1)/R1^2

disp('R1=R2=0')
FA^2*FB*(3*R12+2*FA*R1-FA*R12+6)/12

% ----- S3_case3_int ----

disp('S3_case3_int')
(R2*expl(2,R1*FB)-R1*expl(2,R2*FB))*expl(1,FA*R1)/(R1^2*R2*(R1-R2))

disp('R1=R12')
(-expl(2,R1*FB)+FB*R1*expl(1,R1*FB))*expl(1,FA*R1)/(R1^3)

disp('R1=0')
FA*expl(2,R2*FB)/(R2^2)

disp('R2=0')
(-expl(3,R1*FB)-expl(3,FA*R1)+expl(3,R1)-FB*R1*expl(2,FA*R1))/(R1^3)

disp('0=R1=R12')
FA*FB^2*(2*R1+2*R2+FA*R1-2*FA*R2+6)/12















