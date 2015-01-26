data=csvread('2014-05-28-copied_from_2014-05-15_BAC_for_import_to_matlab.csv');
subplot(2,2,1)
hist(data(:,1),30)
title('Sfp only')
subplot(2,2,3)
hist(data(:,2),30)
title('AcpS only')
print -dpdf '2014-05-28-plot-of-old-normalization.pdf'
