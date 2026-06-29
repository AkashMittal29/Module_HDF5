
% clear;

info=h5info('field.h5');

iter1 = h5read('field.h5','/iter1');
iter2 = h5read('field.h5','/iter2');

value1 = h5read('field.h5','/value1');
value2 = h5read('field.h5','/value2');


%%
iter1_old = iter1;
iter2_old = iter2;
value1_old = value1;
value2_old = value2;


%%
plotH5Hierarchy('field.h5');
