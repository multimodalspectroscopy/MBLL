%plot concentrations

figure
plot(time,conc(:,1),'r')
hold on
plot(time,conc(:,2),'b')
plot(time,conc(:,3),'g')
xlabel('Time (s)')
ylabel('Change in concentration (mM)')
legend('HbO_2','HHb','oxCCO')