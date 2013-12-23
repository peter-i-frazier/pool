a=csvread('forplot.csv',1,1);
h=figure;
hold;
plot(1:100,a(1,:),'r');
plot(1:100,a(2,:),'g')
plot(1:100,a(3,:),'b')
xlabel('No. of recommendations','FontSize',30);
ylabel('P(shortest hit <= 12)','FontSize',30);
axis([0 100 0 1.2]);
legend('Probability of improvement algorithm','Pick most probable from posterior','Mutate from known hits');
print(h,'-dpdf','maxP_comparison.pdf');