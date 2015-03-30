clear all;
clc;

load('pointTargetData.mat');

data = veraStrct.data(80:end,:,:);
fs = 20e6;
[rows_data,col_data, z_data] = size(data);
time = [0:1:rows_data-1]*(1/fs);
time_upsample = [0:1/4:rows_data-1]*(1/fs);

%part a
figure;
imagesc(data(:,:),[min(min(min(data))) max(max(max(data)))]);
colormap('gray');
title('Channel data (pointTargetData.mat)');
%zoom

for kk = 1:128;
        data_upsample(:,:,kk) = interp1(time,data(:,:,kk),time_upsample,'linear');
end

fs_upsample = 4*fs;
speed = 1540; %m/s in body
pixel_size_through_depth = speed/fs_upsample; 
time_zero = 80;

channel = [[-63.5:1:63.5]];

depth = 0.04; %m    
for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    d(ii) = (xe(ii)^2+depth^2)^0.5;
    time_to_point(ii) = d(ii)/speed;
end
time_from_zero = time_to_point(64);
time_from_zero_v = ones(1,length(time_to_point))*time_from_zero;
time_delay = time_to_point - time_from_zero_v;
pixels_delay = time_delay*fs_upsample;
pixels_delay_rounded = round(pixels_delay);

[rows_data col_data z_data] = size(data_upsample);
data_shifted = zeros(rows_data, col_data, z_data);
for b = 1:128
    for a = 1:128
        length_remaining = rows_data-pixels_delay_rounded(a);
        data_shifted(1:length_remaining,a,b) = data_upsample((pixels_delay_rounded(a)+1:end),a,b);
    end
end
% for b = 1:128
%     for a = 1:128
%         length_remaining = rows_data-pixels_delay_rounded(a);
%         data_shifted(pixels_delay_rounded(a)+pixels_delay_rounded(a):end,a,b) = data_upsample((pixels_delay_rounded(a):pixels_delay_rounded(a)+9528),a,b);
%     end
% end

%part b
figure;
imagesc(data_upsample([3500:4500],:,64));
colormap('gray');
figure;
imagesc(data_shifted([3500:4500],:,64));
colormap('gray');
title('Point target, channel 64 (pointTargetData.mat)');
