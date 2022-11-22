%% Jiaheng Yu, BME college, Beihang University
%% data import
d=importdata("16channel recordings.dat");
d1=d.data;
data=d1(:,1); 
test_data=data(1000:10999);%data import, select 10000 data points as sample
%% filter design and testing
filtered_data=my_filter(test_data);
figure(1)
subplot(411)
plot(test_data);
xticklabels (0:0.05:0.5)
title('Raw Data'),xlabel('time(s)')
subplot(412)
plot(filtered_data);
xticklabels (0:0.05:0.5)
title('Filtered Data'),xlabel('time(s)')%raw data and filtered signal
subplot(413)
TEST_DATA=fft(test_data);
stem(abs(TEST_DATA(1:100:10000)));
axis([0 (200000/10000) 0 max(abs(TEST_DATA))]);%spectrum analysis
xlabel('Frequency(Hz)'),ylabel('Amlplitude(dB)')
title('Frequency of Original Data')
subplot(414)
F_DATA=fft(filtered_data);
stem(abs(F_DATA(1:100:10000)));
axis([0 (200000/10000) 0 max(abs(F_DATA))]);
xlabel('Frequency(Hz)'),ylabel('Amlplitude(dB)')
title('Frequency of Filtered Data')
%% data filter
save_data=zeros(160000,16);
for cnt=1:16
    data=d1(:,cnt);
    test_data=data;
    filtered_data=my_filter(test_data);
    save_data(:,cnt)=filtered_data;
end
figure (2)
for cnt=1:16
    subplot(16,1,cnt)
    plot(save_data(:,cnt));
    %axis off
end
%% time latency
save_max=zeros(1,16);
for cnt=1:16
    data=d1(:,cnt);
    test_data=data(15000:18999);%data import, select 10000 data points as sample
    [m,post]=max(test_data);
    save_max(cnt)=post;
end
pos_c=[16 10 5 3 8 1 2 6 12 11 15 13 9 7 14 4];
x=[4 4 3 4 3 4 3 3 2 2 1 2 2 1 1 1];
y=[1 2 1 3 2 4 3 4 1 2 1 3 4 2 3 4];%map the device position with channel number
time=zeros(4,4);
min_t=min(save_max);
for cnt=1:16
    time(y(cnt),x(cnt))=(save_max(pos_c(cnt))-min_t)/20;
end
%time=rot90(time,3);
x=1:4;y=1:4;
xx=1:0.1:4;yy=1:0.1:4;
[X,Y]=meshgrid(xx,yy);
time_new=interp2(x,y,time,X,Y,'cubic');%interpration
figure(3)
contourf(time_new(1:30,1:30),10);%isochronal map of time latency
axis equal
xticks(0:10:30)
xticklabels([1 2 3 4])
yticks(0:10:30)
yticklabels([1 2 3 4])
importdata('CBR_wet.rgb');
grid
colorbar
color = ncl_colormap('CBR_wet');
colormap(color)
colorbar
%% filter function
function y = my_filter(x)

% MATLAB Code
% Generated by MATLAB(R) 9.10 and DSP System Toolbox 9.12.
% Generated on: 18-Nov-2022 23:40:49

%#codegen

persistent Hd;

if isempty(Hd)
    
    % filter design codes:
    % % Equiripple Lowpass filter designed using the FIRPM function.
    %
    % % All frequency values are in Hz.
    % Fs = 20000;  % Sampling Frequency
    %
    % N     = 10;    % Order
    % Fpass = 0;     % Passband Frequency
    % Fstop = 6000;  % Stopband Frequency
    % Wpass = 1;     % Passband Weight
    % Wstop = 80;    % Stopband Weight
    % dens  = 20;    % Density Factor
    %
    % % Calculate the coefficients using the FIRPM function.
    % b  = firpm(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop], ...
    %            {dens});
    
    Hd = dsp.FIRFilter( ...
        'Numerator', [0.00265709699723654 0.0173960054136978 ...
        0.0572568013165367 0.123007522938122 0.18957134876703 0.218137374019024 ...
        0.18957134876703 0.123007522938122 0.0572568013165367 0.0173960054136978 ...
        0.00265709699723654]);
end

y = step(Hd,double(x));
end
%% colormap functions
function color = ncl_colormap(colorname)

temp = import_ascii([colorname '.rgb']);
temp(1:2) = [];
temp = split(temp,'#');
temp = temp(:,1);
% color = deblank(color);
temp = strtrim(temp);
temp = regexp(temp, '\s+', 'split');
for i=1:size(temp,1)
    color(i,:) = str2double(temp{i});    
end
color = color/255;
end
% Edited Time:2019-02-22
function ascii = import_ascii(file_name)
i = 1;
fid = fopen(file_name);
while feof(fid) ~= 1
    tline = fgetl(fid);
    ascii{i,1} = tline; i = i + 1;
end
fclose(fid);
end