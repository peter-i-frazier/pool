a=csvread('PP_use_newData#2.csv',1,1);
h=figure;
hold;
plot(1:100,a(1,:),'r');
plot(1:100,a(2,:),'g')
plot(1:100,a(3,:),'b')
xlabel('No. of recommendations','FontSize',30);
ylabel('P(shortest hit <= 12)','FontSize',30);
axis([0 100 0 1.2]);
legend('maxP','Ranking estimated probability of hit','mutate');
print(h,'-dpng','maxP_comparison.png');