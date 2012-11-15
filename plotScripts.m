% figure
% hist(vDP{analyze(1)},20)
% hold on
% hist(-vSmooth{analyze(1)},20)
% pause
% close all
% 
% figure
% subplot(1,2,1)
% rose(theta{analyze(1)})
% legend('Original theta')
% subplot(1,2,2)
% rose(thetaDP{analyze(1)})
% legend('DP theta')
% pause
% close all
% 
% clear vDPSegX vDPSegY


%% Start plots

% %% Plot trajectories, color coded by transport state
% figure
% for ii = analyze
%     plotBySeg(xPos{ii}, yPos{ii},event{ii},segState{ii})
%     title('Plot by segment transport state')
%     pause
%     %plotByDir(xPos{ii},yPos{ii},direction{ii})
%     %pause
%     clf
% end

%% Plot by mean segment direction
% for ii = analyze
%     plotBySegDir(xPos{ii},yPos{ii},event{ii},segDir{ii})
%     title('Plot by segment direction, green away from body')
%     pause
%     clf
% end

% 
% hist(meanRunV(analyze)./meanStagV(analyze))
% pause

%% Plot histogram of distance travelled in active segments
% for ii = analyze
%     hist(segDistance{ii}(find(segState{ii}==3)))
%     pause
% end
clf
%% Plot 'bimodal' histogram of active velocities, one for direction=1, other for direction=0:
% for ii = analyze
%     hist(vSmooth{ii}(find(state{ii}==3 & direction{ii}==1)))
%     hold on
%     hist(-vSmooth{ii}(find(state{ii}==3 & direction{ii}==0)))
%     pause
% end

% for ii = analyze
%     Q = vSmooth{ii}(find(state{ii}==3 & direction{ii}==1));
%     R = vSmooth{ii}(find(state{ii}==3 & direction{ii}==0));
%     histPair(Q,R)
%     pause
%     clf
% end

close all