w=20;

load(strcat('dosNV200Npsi200w',num2str(w),'.mat'));
save(strcat('int1w',num2str(w),'.dat'),'int1','-ascii');
save(strcat('int2w',num2str(w),'.dat'),'int2','-ascii');
peakdiff;
save(strcat('peakdiffw',num2str(w),'.dat'),'peakd','-ascii');
