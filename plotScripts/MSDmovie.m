%% MSD Movie
for ii = 300:1458
    loglog(tau{2}{ii},MSD{2}{ii},'b',tau{6}{ii},MSD{6}{ii},'r')
    axis([0 1 1E-16 1E-14])
    pause
end