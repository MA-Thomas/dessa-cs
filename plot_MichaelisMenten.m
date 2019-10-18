% PLOT RESULTS - Michaelis-Menten

indices = timePoints_trials{1,1} > 0;
times = timePoints_trials{1,1}(indices);
E = numE_trials{1,1}(indices);
S = numS_trials{1,1}(indices);
ES = numES_trials{1,1}(indices);
P = numP_trials{1,1}(indices);

figure
plot(times,E)
hold on
plot(times,S)
plot(times,ES)
plot(times,P)
xlabel('Time (s)')
ylabel('Molecules')
legend('E','S','ES','P')
ylim([0 1000])
xlim([0 100])
hold off
% title(['THARD-preCompute,THARD-initializ,THARD-main: ',num2str(8),', ',num2str(8),', ',num2str(8),...
%     ' Code to restrict partners uncommented'])