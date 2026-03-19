
% Modify the following commands as per the hierarchy of the h5 file 
% as defined in the Fortran main code.
clear;
info=h5info('field.h5');

% x = h5read('field.h5','/subgrid/x');
% y = h5read('field.h5','/subgrid/y');

data = h5read('field.h5','/iter');
x = h5read('field.h5','/grid/x');
y = h5read('field.h5','/grid/y');
vel = h5read('field.h5','/field/velocity');
pres = h5read('field.h5','/field/pressure');

points = h5read('marker.h5','/points');


%% Reading a hyperslab
% Modify the following commands as per info.
start = [1,1,1]; % Starting indices in each dimension.
count = info.Groups.Datasets.ChunkSize; % No. of elements to read in each dimension.
data2 = h5read('velocity.h5','/mydata/velocity',start, count);

