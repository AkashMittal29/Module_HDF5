
info=h5info('field.h5');

x = h5read('field.h5','/grid/x');
y = h5read('field.h5','/grid/y');
data1 = h5read('field.h5','/iter');
vel = h5read('field.h5','/field/velocity');
pres = h5read('field.h5','/field/pressure');
points = h5read('marker.h5','/points');

% Reading a hyperslab
start = [1,1,1]; % Starting indices in each dimension.
count = info.Groups.Datasets.ChunkSize; % No. of elements to read in each dimension.
data2 = h5read('velocity.h5','/mydata/velocity',start, count);

