clc

FA=0.5;
npts=100;
pts=logspace(-8,2,npts)*(-1);
data=zeros(npts,1);
Old=zeros(npts,1);

for ii=1:npts
    R1=0.3;
    R2=pts(ii);
    data(ii)=S3_case3_int(FA,R1,R2);
    %Old(ii)=S3_case3_intOld(FA,R1,R2);
end
figure
semilogx(pts,data,pts,Old,'+')
ylabel('S3_case2_int')
xlabel('R2')
title('R1=0.3, FA=0.5')

for ii=1:npts
    R1=0.000003;
    R2=pts(ii);
    data(ii)=S3_case3_int(FA,R1,R2);
    %Old(ii)=S3_case3_intOld(FA,R1,R2);
end
figure
semilogx(pts,data,pts,Old,'+')
ylabel('S3_case2_int')
xlabel('R2')
title('R1=0.000003, FA=0.5')

for ii=1:npts
    R1=0.3;
    R2=R1+pts(ii);
    data(ii)=S3_case3_int(FA,R1,R2);
    %Old(ii)=S3_case3_intOld(FA,R1,R2);
end
figure
semilogx(pts,data,pts,Old,'+')
ylabel('S3_case2_int')
xlabel('R2-R1')
title('R1=0.3, FA=0.5')


% for ii=1:npts
%     R1=0.3;
%     R12=pts(ii);
%     R3=0.301;
%     data(ii)=S4_case3_int(FA,R1,R12,R3);
%     Old(ii)=S4_case3_intOld(FA,R1,R12,R3);
% end
% figure
% semilogx(pts,data,pts,Old,'+')
% ylabel('S3_case2_int')
% xlabel('R2')
% title('R1=0.3, FA=0.5')
% 
% for ii=1:npts
%     R1=0.3;
%     R12=pts(ii);
%     R3=0.4;
%     data(ii)=S4_case3_int(FA,R1,R12,R3);
%     Old(ii)=S4_case3_intOld(FA,R1,R12,R3);
% end
% figure
% semilogx(pts,data,pts,Old,'+')
% ylabel('S3_case2_int')
% xlabel('R1-R12')
% title('R1=0.3, FA=0.5')
% ylim([0,1])