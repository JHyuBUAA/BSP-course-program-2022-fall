%% BSP course project, Jiaheng Yu, BME college, Beihang University, Beijing.
%% part1
%% data import
d=importdata("16channel recordings.dat");
d1=d.data;
data=d1(:,1); 
test_data=data(1000:10999);%data import, select 10000 data points as sample
test_data=test_data-mean(test_data(1000:3000));
%% filter design and testing
filtered_data=my_filter(test_data);
figure(1)
subplot(411)
plot(test_data);
xticklabels (0:0.05:0.5)
title('Raw Data'),xlabel('time(s)')
l=length(test_data);
f=2*pi/l:2*pi/l:2*pi;
subplot(412)
plot(filtered_data);
xticklabels (0:0.05:0.5)
title('Filtered Data'),xlabel('time(s)')%raw data and filtered signal
subplot(413)
TEST_DATA=fft(test_data);
plot(f,abs(TEST_DATA));
axis([0 (200000/10000) 0 max(abs(TEST_DATA))]);%spectrum analysis
xlabel('Frequency(Hz)'),ylabel('Amlplitude(dB)')
title('Frequency of Original Data')
subplot(414)
F_DATA=fft(filtered_data);
plot(f,abs(F_DATA));
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
%% part2
%% parameter definition
% import data
data=reshape(d1(:,2:6),1,160000*5);
test_data =data;

thres = 5; %define multiple of sigma (for threshold settings)
threshold = thres.*median(abs(test_data)/0.6745);%Obey the standard normal distribution, so that the error probability falls in the 25% to 75% area

sf = 20; % define sample frequency  20000
clusternum = 3; %define number of cluster in analysis
pre_time =3; post_time = 3; %in ms, acquisition time before (pre_time) and after (post_time)detection of a waveform peak
time_stamp = [];
waveform = [];
ii = pre_time*20;

%% spike detection
count = 0;
length_data=length(test_data);
while ii < length_data
    tmp = test_data(ii);
    if tmp > threshold
        if post_time*20+ii < length(test_data)
            count = count + 1;
            time_stamp(count) = ii;
            waveform(count,:) = test_data((-pre_time*20:post_time*20)+ii);%save the datas pre_time*30 before and post_time*30 after the detected spike, 30 data points are recorded in a ms
        end
        ii = ii + post_time * 20;
    else
        ii = ii + 1;
    end
end
%disp(time_stamp);
[row,column] = size(waveform);
[row2,column2] = size(data);%number of the detected waveform
figure
subplot(2,1,1)
plot((-pre_time*20:post_time*20)/20,waveform');
title(['total spike number is ',num2str(row)]);
xlabel('time(ms)');
ylabel('Voltage(uV)');

subplot(2,1,2)
plot((-pre_time*20:post_time*20)/20,mean(waveform),LineWidth=2);
title('mean waveform');
xlabel('time(ms)');
ylabel('Voltage(uV)');

figure
time_for_plot = (1:length(test_data))/sf;% sf: sampling frequency
subplot(2,1,1)
plot(time_for_plot,test_data);
xlim([0,length(data)/sf]);
peaks=findpeaks(test_data,'MinPeakHeight',-threshold);
title('Filtered data');
xlabel('time(ms)');
ylabel('Voltage(uV)');

subplot(2,1,2)
hold on
for ii = 1:length(time_stamp)
plot([time_stamp(ii),time_stamp(ii)]/sf,[-1,1],'k');
end
xlim([0,length(data)/sf]);
ylim([-10,10]);
title('Raster plot');
xlabel('time(ms)');
ylabel('Raster');

%% dimension reduction (PCA)
%waveform = zscore(waveform); % Standardized data
figure
[coeff,score,latent,tsquared,explained] = pca(waveform);
h = biplot(coeff(:,1:2),'Scores',score(:,1:2));
xlabel('First PC');
ylabel('Second PC');
title('Principal Component Analysis');

% pca plot
% mapcaplot(waveform);

%% cluster (kmean)
figure
bar(explained)
title('Explained Variance')
ylabel('PC')

% Retain first two principal components
yeastPC = score(:,1:2);
figure
[clusters, centroid] = kmeans(yeastPC,clusternum);% clusternum=2
gscatter(yeastPC(:,2),yeastPC(:,1),clusters)
xlabel('First PC');
ylabel('Second PC');
title('Principal Component Scatter Plot with Colored Clusters');

%% Plot result
figure
col=[0.86,0.43,0.34;0.43,0.68,0.82;0.72,0.13,0.19;0.06,0.27,0.5;];
for c = 1:clusternum
    subplot(3,clusternum,c);
    plot((-pre_time*20:post_time*20)/20,waveform((clusters == c),:)','Color',col(c,:));
    xlabel('time (ms)');
    ylabel('Voltage (uV)');
    title(['Cluster',num2str(c)]);
    
    subplot(3,clusternum,c+clusternum);
    plot((-pre_time*20:post_time*20)/20,mean(waveform((clusters == c),:)),'Color',col(c,:),'Linewidth',2);
    title('Mean spikes');
   
    x=(-pre_time*20:post_time*20)/20;
    curve1=mean(waveform((clusters == c),:))-std(waveform((clusters == c),:));
    curve2=mean(waveform((clusters == c),:))+std(waveform((clusters == c),:));
    subplot(3,clusternum,c+2*clusternum);
    plot(x, curve1, 'k--', 'LineWidth', 1,'Color',col(c,:));
    hold on;
    plot(x, curve2, 'k--', 'LineWidth', 1,'Color',col(c,:));
    hold on;                                    % add Paul
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [0.85 0.85 0.85]);
    axis tight   
    title('Spike group (mean+-std)');
    % add mean spike
    plot((-pre_time*20:post_time*20)/20,mean(waveform((clusters == c),:)),'Linewidth',2,'Color',col(c,:));
    
    %min dVdt
    min(diff(mean(waveform((clusters == c),:)))*20);
end
% suptitle('Clustering of Profiles');

%% Color raster plot
figure
for i=1:clusternum
    Snumb(i)=0;%counting spikes in each cluster
end

for ii = 1:length(clusters)
    if clusters(ii)==1
        Snumb(1)=Snumb(1)+1;
        time_stamps(1,Snumb(1))=time_stamp(ii);
        plot([time_stamp(ii),time_stamp(ii)]/sf,[-1,1],'Linewidth',2,'Color',col(1,:));
        hold on
    elseif clusters(ii)==2
        Snumb(2)=Snumb(2)+1;
        time_stamps(2,Snumb(2))=time_stamp(ii);
        plot([time_stamp(ii),time_stamp(ii)]/sf,[-1,1],'Linewidth',2,'Color',col(2,:));
        hold on
    elseif clusters(ii)==3
     Snumb(3)=Snumb(3)+1;
        time_stamps(3,Snumb(3))=time_stamp(ii);
        plot([time_stamp(ii),time_stamp(ii)]/sf,[-1,1],'Linewidth',2,'Color',col(3,:));
        hold on
    elseif clusters(ii)==4
        Snumb(4)=Snumb(4)+1;
        time_stamps(4,Snumb(4))=time_stamp(ii);
        plot([time_stamp(ii),time_stamp(ii)]/sf,[-1,1],'Linewidth',2,'Color',col(4,:));
        hold on
    elseif clusters(ii)==5
        Snumb(5)=Snumb(5)+1;
        time_stamps(5,Snumb(5))=time_stamp(ii);
        plot([time_stamps(ii),time_stamp(ii)]/sf,[-1,1],'Linewidth',2,'Color',col(5,:));
        hold on
        end
end
xlim([0,length(data)/sf]);
ylim([-10,10]);
title('Colored Raster plot');
xlabel('time(s)');
ylabel('Raster');

%% Additional analysis for each cluster
figure()
for c = 1:clusternum
    subplot(clusternum,1,c);
        %Half-width
        [pks,locs,widths,proms]=findpeaks(mean(waveform((clusters == c),:)),(-pre_time*20:post_time*20)/20,'Annotate','extents','WidthReference','halfheight');
        %disp(num2str(pks));
        %disp(num2str(widths));
        title('Signal Peak Widths')
        hold on
        curve1=mean(waveform((clusters == c),:))-std(waveform((clusters == c),:));
        curve2=mean(waveform((clusters == c),:))+std(waveform((clusters == c),:));
        findpeaks(curve1,(-pre_time*20:post_time*20)/20,'Annotate','extents','WidthReference','halfheight');
        hold on
        findpeaks(curve2,(-pre_time*20:post_time*20)/20,'Annotate','extents','WidthReference','halfheight');
        %legend('off');
        hold on
        %ISI calculation
            for j=2:Snumb(c)
                if time_stamps(c,j)>0
                    ISI(c,j)=time_stamps(c,j)-time_stamps(c,j-1);
                end  
            end
end
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
function ascii = import_ascii(file_name)
i = 1;
fid = fopen(file_name);
while feof(fid) ~= 1
    tline = fgetl(fid);
    ascii{i,1} = tline; i = i + 1;
end
fclose(fid);
end
%% present with love to Shushu Chen