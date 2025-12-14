clc
clear all
close all

filename = 'D:\Chrome Download\2.nc';
info = ncinfo(filename);
disp(info)

for i = 1:length(info.Variables)
    disp(info.Variables(i).Name)
end

%%
% time = ncread(filename,'time');
% t0 = datetime(1970,1,1,'TimeZone','UTC'); % 基准时间
% time_dt = t0 + seconds(time);
% 
% % 打印前 10 行
% fprintf('     Original Seconds          Converted Time (UTC)\n');
% fprintf('-------------------------------------------------------------\n');
% for i = 1:10
%     fprintf('%15.0f     %s\n', time(i), datestr(time_dt(i)));
% end

%%
% 读取变量
time = ncread(filename,'time');
VMDR_SW1  = ncread(filename,'VMDR_SW1');


% 正确基准时间：Unix epoch
t0 = datetime(1970,1,1,'TimeZone','UTC');
time_dt = t0 + seconds(time);


fprintf('Time (UTC)            主涌浪方向SW1(deg)      \n');
fprintf('-------------------------------------------------------------------------------\n');
for i = 1:233 % 969
    fprintf('%s        %8.2f   \n', ...
        datestr(time_dt(i)), VMDR_SW1(i));
end




