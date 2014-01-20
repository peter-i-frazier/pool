% This datafile was originally called PP_use_newData#2.csv
% It contains P(shortest hist <= 12) vs. no. of recommendations for three
% different methods,
% maxP is choosing the next thing to add to maximize probability of
% improvement, adding 10 at a time;
% ranking estimated probability;
% mutation, averaged over 100 independent replications, where in each
% replication S contains the first |S| of a set of 100 randomly mutated
% peptides.
a=csvread('forplot.csv',1,1); 
h=figure;
hold;
plot(1:100,a(1,:),'r');
plot(1:100,a(2,:),'g')
plot(1:100,a(3,:),'b')
xlabel('No. of recommendations','FontSize',30);
ylabel('P(shortest hit <= 12)','FontSize',30);
axis([0 100 0 1.2]);
legend('maxP','Ranking estimated probability of hit','mutate');
% Or, if you have room for a more verbose legend:
% legend('Probability of improvement algorithm','Pick most probable from posterior','Mutate from known hits');
print(h,'-dpng','maxP_comparison.png');
