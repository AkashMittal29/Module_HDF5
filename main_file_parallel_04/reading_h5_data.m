
clear;

info=h5info('field.h5');

x_full = h5read('field_full.h5','/grid/x');
y_full = h5read('field_full.h5','/grid/y');

iter = h5read('field.h5','/iter');
x = h5read('field.h5','/subgrid/x');
y = h5read('field.h5','/subgrid/y');
vel = h5read('field.h5','/field/velocity');
pres = h5read('field.h5','/field/pressure');


%% Checks
i = 5;
vel_actual = 0.01*ones(size(x(:,:)))*double(iter(i))+x(:,:).^2+y(:,:);
sum(sum(vel_actual-vel(:,:,iter(i))))

pres_actual = 0.01*ones(size(x(:,:)))*double(iter(i))+x(:,:).^3+y(:,:).^2;
sum(sum(pres_actual-pres(:,:,iter(i))))


%% Reading a hyperslab
start = [1,1,1]; % Starting indices in each dimension.
count = info.Groups(1).Datasets.ChunkSize; % No. of elements to read in each dimension.
data2 = h5read('field.h5','/field/velocity',start, count);

