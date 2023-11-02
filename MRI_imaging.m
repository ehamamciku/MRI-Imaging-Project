%% Loading the data files
clc; clear all; close all;

% selecting the data files to be used
% first select the good data files then the bad data files
[file1,path1] = uigetfile('*.mat');
[file2,path2] = uigetfile('*.mat');
[file3,path3] = uigetfile('*.mat');
[file4,path4] = uigetfile('*.mat');
[file5,path5] = uigetfile('*.mat');
[file6,path6] = uigetfile('*.mat');

good_data1 = cell2mat(struct2cell(load(fullfile(path1, file1))));
good_data2 = cell2mat(struct2cell(load(fullfile(path2, file2))));
good_data3 = cell2mat(struct2cell(load(fullfile(path3, file3))));
bad_data1 = cell2mat(struct2cell(load(fullfile(path4, file4))));
bad_data2 = cell2mat(struct2cell(load(fullfile(path5, file5))));
bad_data3 = cell2mat(struct2cell(load(fullfile(path6, file6))));

%change the number 49 depending on the structure of the directory
slice_number = path1(49);

%% Detection of the outliers

dif1 = good_data1 - bad_data1;
dif2 = good_data2 - bad_data2;
dif3 = good_data3 - bad_data3;

magdif1 = abs(dif1);
magdif2 = abs(dif2);
magdif3 = abs(dif3);

for i = 1:size(magdif1,2)
    pow_col1(i) = sqrt(sum(magdif1(:,i).*magdif1(:,i)))/512;  
    pow_col2(i) = sqrt(sum(magdif2(:,i).*magdif2(:,i)))/512;  
    pow_col3(i) = sqrt(sum(magdif3(:,i).*magdif3(:,i)))/512;  
end

[TF1,P1] = islocalmax(pow_col1);
P1(TF1);
[TF2,P2] = islocalmax(pow_col2);
P2(TF2);
[TF3,P3] = islocalmax(pow_col3);
P3(TF3);

threshold1 = P1(6);
threshold2 = P2(6);
threshold3 = P3(6);

I1 = find(P1>=threshold1);
I2 = find(P2>=threshold2);
I3 = find(P3>=threshold3);

Inew1 = [I1(1:6),I1(9:end)];
Inew2 = [I2(1:7),I2(10:end)];
Inew3 = [I3(1:3),I3(5:end)];

outlier_data1 = bad_data1;

for i = Inew1
    for j = 1:512
        outlier_data1(j,i) = 255;
    end
end

%% Linear Interpolation in k-space

linear_restored1 = bad_data1;
linear_restored2 = bad_data2;
linear_restored3 = bad_data3;

% Interpolating seperately for each channel data

for i = 1:length(Inew1)
    col1 = Inew1(i);
    linear_restored1(:,col1) = (linear_restored1(:,col1-1)+linear_restored1(:,col1+1))/2;
end

for i = 1:length(Inew2)
    col2 = Inew2(i);
    linear_restored2(:,col2) = (linear_restored2(:,col2-1)+linear_restored2(:,col2+1))/2;
end

for i = 1:length(Inew3)
    col3 = Inew3(i);
    linear_restored3(:,col3) = (linear_restored3(:,col3-1)+linear_restored3(:,col3+1))/2;
end

%% Cubic Interpolation in k-space

cubic_restored1 = bad_data1;
cubic_restored2 = bad_data2;
cubic_restored3 = bad_data3;

cubic_restored1(:,Inew1) = 0;
cubic_restored2(:,Inew2) = 0;
cubic_restored3(:,Inew3) = 0;

cubic_filt = [-0.0625 0 0.5625 1 0.5625 0 -0.0625];

% Interpolating seperately for each channel data

for i = 1:length(Inew1)
    for j = 1:size(cubic_restored1,1)
        coll1 = Inew1(i);
        seg1 = conv(cubic_filt,cubic_restored1(j,coll1-5:coll1+5),'same');
        cubic_restored1(j,coll1) = seg1(4);
    end
end

for i = 1:length(Inew2)
    for j = 1:size(cubic_restored2,1)
        coll2 = Inew2(i);
        seg2 = conv(cubic_filt,cubic_restored2(j,coll2-5:coll2+5),'same');
        cubic_restored2(j,coll2) = seg2(4);
    end
end

for i = 1:length(Inew3)
    for j = 1:size(cubic_restored3,1)
        coll3 = Inew3(i);
        seg3 = conv(cubic_filt,cubic_restored3(j,coll3-5:coll3+5),'same');
        cubic_restored3(j,coll3) = seg3(4);
    end
end

%% Using a Gaussian filter in k-space

gauss_restored1 = bad_data1;
gauss_restored2 = bad_data2;
gauss_restored3 = bad_data3;

% Defining the Gaussian filter
filt_len = 3;
sigma = 0.8;
n = (filt_len-1)/2;
for i = -n:1:n
    for j = -n:1:n
    gaus_filter(i+n+1,j+n+1) = exp(((i)^2 + (j)^2)/(2*sigma*sigma));
    end
end

% Normalizing the impulse response of the filter
K = sum(gaus_filter,'all');
gaus_filter = gaus_filter/K;

% Applying the Gaussian filter

for i = 1:length(Inew1)
    colll1 = Inew1(i);
    segg1 = conv2(gauss_restored1(:,colll1-1:colll1+1),gaus_filter,'same');
    gauss_restored1(:,colll1) = segg1(:,2);
end

for i = 1:length(Inew2)
    colll2 = Inew2(i);
    segg2 = conv2(gauss_restored2(:,colll2-1:colll2+1),gaus_filter,'same');
    gauss_restored2(:,colll2) = segg2(:,2);
end
for i = 1:length(Inew3)
    colll2 = Inew3(i);
    segg3 = conv2(gauss_restored3(:,colll2-1:colll2+1),gaus_filter,'same');
    gauss_restored3(:,colll2) = segg3(:,2);
end

%% Using a Low-Pass filter in k-space

lowpass_restored1 = bad_data1;
lowpass_restored2 = bad_data2;
lowpass_restored3 = bad_data3;

% Defining the Low-Pass filter
filt_len = 3;
n = (filt_len-1)/2;
for i = -n:1:n
    for j = -n:1:n
    lowpass_filter(i+n+1,j+n+1) = 1;
    end
end

% Normalizing the impulse response of the filter
K = sum(lowpass_filter,'all');
lowpass_filter = lowpass_filter/K;

% Applying the Low-Pass filter

for i = 1:length(Inew1)
    colll1 = Inew1(i);
    segg1 = conv2(gauss_restored1(:,colll1-1:colll1+1),lowpass_filter,'same');
    lowpass_restored1(:,colll1) = segg1(:,2);
end

for i = 1:length(Inew2)
    colll2 = Inew2(i);
    segg2 = conv2(gauss_restored2(:,colll2-1:colll2+1),lowpass_filter,'same');
    lowpass_restored2(:,colll2) = segg2(:,2);
end
for i = 1:length(Inew3)
    colll2 = Inew3(i);
    segg3 = conv2(gauss_restored3(:,colll2-1:colll2+1),lowpass_filter,'same');
    lowpass_restored3(:,colll2) = segg3(:,2);
end

%% Using a Wiener filter in k-space

% Filter initializations
filterorder = 3;
nn = (filterorder-1)/2;

% Autocorrelation and cross-correlation will be calculated upon spatial eye images
% Translating the k-space images into spatial eye images

% 1. X - dimension of the K-Space data    - 128
% 2. Y - dimension of the K-Space data    - 512

% IFFT of k-space data for good data
%channel 1 
Data_img_good_wiener(:,:,1) = ifftshift(ifft2(good_data1),1);
%channel 2
Data_img_good_wiener(:,:,2) = ifftshift(ifft2(good_data2),1);
%channel 3
Data_img_good_wiener(:,:,3) = ifftshift(ifft2(good_data3),1);

% IFFT of k-space data for bad data
%channel 1 
Data_img_bad_wiener(:,:,1) = ifftshift(ifft2(bad_data1),1);
%channel 2
Data_img_bad_wiener(:,:,2) = ifftshift(ifft2(bad_data2),1);
%channel 3
Data_img_bad_wiener(:,:,3) = ifftshift(ifft2(bad_data3),1);

% Matrix initializations
M = size(Data_img_bad_wiener, 2);
N = size(Data_img_bad_wiener, 1) - 1;

for j = 1:size(Data_img_bad_wiener, 3)
    % Autocorrelation of the bad data
    rx = zeros(filterorder, size(Data_img_bad_wiener, 2));
    for k = 0:filterorder - 1
        for n = k:N
            rx(k+1, :) = rx(k+1, :) + Data_img_bad_wiener(n+1,:,j) .* conj(Data_img_bad_wiener(n-k+1,:,j));
        end
    end
    % Crosscorrelation between bad data and good data
    rdx = zeros(filterorder, size(Data_img_bad_wiener, 2));
    for k = 0:filterorder - 1
        for n = k:N
            rdx(k+1, :) = rdx(k+1, :) + Data_img_good_wiener(n+1,:,j) .* conj(Data_img_bad_wiener(n-k+1,:,j));
        end
    end
    % Wiener-Hopf equations
    w = zeros(size(rx));
    for i = 1:M
        Rx = toeplitz(rx(:, i)).';
        w(:, i) = Rx \ rdx(:, i);
    end
    % Normalization of the filter coefficients
    w(:, i) = w(:, i) / max(abs(w(:, i)));
    wiener_restored1 = bad_data1;
    wiener_restored2 = bad_data2;
    wiener_restored3 = bad_data3;
    filter1 = w(1,:);
    filter2 = w(2,:);
    filter3 = w(3,:);
    for i = 1:length(Inew1)
        colll1 = Inew1(i);
        segg1 = conv2(bad_data1(:,colll1-nn:colll1+nn),filter1,'same');
        wiener_restored1(:,colll1) = segg1(:,2);
    end
    for i = 1:length(Inew2)
        colll2 = Inew2(i);
        segg2 = conv2(bad_data2(:,colll2-nn:colll2+nn),filter2,'same');
        wiener_restored2(:,colll2) = segg2(:,2);
    end
    for i = 1:length(Inew3)
        colll3 = Inew3(i);
        segg3 = conv2(bad_data3(:,colll3-nn:colll3+nn),filter3,'same');
        wiener_restored3(:,colll3) = segg3(:,2);
    end
end

%% Translating the k-space images into spatial eye images

% 1. X - dimension of the K-Space data    - 128
% 2. Y - dimension of the K-Space data    - 512

% IFFT of k-space data for good data
%channel 1 
Data_img_good(:,:,1) = ifftshift(ifft2(good_data1),1);
%channel 2
Data_img_good(:,:,2) = ifftshift(ifft2(good_data2),1);
%channel 3
Data_img_good(:,:,3) = ifftshift(ifft2(good_data3),1);

% IFFT of k-space data for bad data
%channel 1 
Data_img_bad(:,:,1) = ifftshift(ifft2(bad_data1),1);
%channel 2
Data_img_bad(:,:,2) = ifftshift(ifft2(bad_data2),1);
%channel 3
Data_img_bad(:,:,3) = ifftshift(ifft2(bad_data3),1);

% IFFT of k-space data for restored bad data with a Gaussian filter
%channel 1 
Data_img_gauss1 = ifftshift(ifft2(gauss_restored1),1);
%channel 2
Data_img_gauss2 = ifftshift(ifft2(gauss_restored2),1);
%channel 3
Data_img_gauss3 = ifftshift(ifft2(gauss_restored3),1);

% IFFT of k-space data for restored bad data with linear interpolation
%channel 1 
Data_img_linear1 = ifftshift(ifft2(linear_restored1),1);
%channel 2
Data_img_linear2 = ifftshift(ifft2(linear_restored2),1);
%channel 3
Data_img_linear3 = ifftshift(ifft2(linear_restored3),1);

% IFFT of k-space data for restored bad data with cubic interpolation
%channel 1 
Data_img_cubic1 = ifftshift(ifft2(cubic_restored1),1);
%channel 2
Data_img_cubic2 = ifftshift(ifft2(cubic_restored2),1);
%channel 3
Data_img_cubic3 = ifftshift(ifft2(cubic_restored3),1);

% IFFT of k-space data for restored bad data with a Low-Pass filter
%channel 1 
Data_img_lowpass1 = ifftshift(ifft2(lowpass_restored1),1);
%channel 2
Data_img_lowpass2 = ifftshift(ifft2(lowpass_restored2),1);
%channel 3
Data_img_lowpass3 = ifftshift(ifft2(lowpass_restored3),1);

% IFFT of k-space data for restored bad data with a Wiener filter
%channel 1 
Data_img_wiener1 = ifftshift(ifft2(wiener_restored1),1);
%channel 2
Data_img_wiener2 = ifftshift(ifft2(wiener_restored2),1);
%channel 3
Data_img_wiener3 = ifftshift(ifft2(wiener_restored3),1);

% Initializations
clear_comp = linspace(10,0.1,size(Data_img_good,2)).^2; 
clear_matrix = repmat(clear_comp,[size(Data_img_good,1) 1]);
crop_x = [128 + 60 : 348 - 33]; % crop coordinates
std_within = 0.995; 

% For original: good data
eye_raw_good = sqrt(abs(squeeze(Data_img_good(:,:,1))).^2 + ...
           abs(squeeze(Data_img_good(:,:,2))).^2 + ...
           abs(squeeze(Data_img_good(:,:,3))).^2).* clear_matrix;  
eye_raw_good = eye_raw_good(crop_x, :);
eye_visualize_good = reshape(squeeze(eye_raw_good(:,:)),[128 128]);   
[aa_good, val_good] = hist(eye_visualize_good(:),linspace(0,max(...
                                    eye_visualize_good(:)),1000));
    thresh_good = val_good(find(cumsum(aa_good)/sum(aa_good) > std_within,1,'first'));
eye_visualize_good = uint16(eye_visualize_good * 65536 / thresh_good); 

% For original: bad data
eye_raw_bad = sqrt(abs(squeeze(Data_img_bad(:,:,1))).^2 + ...
           abs(squeeze(Data_img_bad(:,:,2))).^2 + ...
           abs(squeeze(Data_img_bad(:,:,3))).^2).* clear_matrix;
eye_raw_bad_before_crop = eye_raw_bad; %for the later section
eye_raw_bad = eye_raw_bad(crop_x, :);
eye_visualize_bad = reshape(squeeze(eye_raw_bad(:,:)),[128 128]);   
[aa_bad, val_bad] = hist(eye_visualize_bad(:),linspace(0,max(...
                                    eye_visualize_bad(:)),1000));
    thresh_bad = val_bad(find(cumsum(aa_bad)/sum(aa_bad) > std_within,1,'first'));
eye_visualize_bad = uint16(eye_visualize_bad * 65536 / thresh_bad);

% For restored: Gaussian filter
eye_raw_gaussian  = sqrt(abs(squeeze(Data_img_gauss1)).^2 + ...
           abs(squeeze(Data_img_gauss2)).^2 + ...
           abs(squeeze(Data_img_gauss3)).^2).* clear_matrix;  
eye_raw_gaussian = eye_raw_gaussian(crop_x, :);
eye_visualize_gaussian = reshape(squeeze(eye_raw_gaussian(:,:)),[128 128]);  
[aa_gauss, val_gauss] = hist(eye_visualize_gaussian(:),linspace(0,max(...
                                    eye_visualize_gaussian(:)),1000));
    thresh_gauss = val_gauss(find(cumsum(aa_gauss)/sum(aa_gauss) > std_within,1,'first'));
eye_visualize_gaussian = uint16(eye_visualize_gaussian * 65536 / thresh_gauss);

% For restored: Low-Pass filter
eye_raw_lowpass  = sqrt(abs(squeeze(Data_img_lowpass1)).^2 + ...
           abs(squeeze(Data_img_lowpass2)).^2 + ...
           abs(squeeze(Data_img_lowpass3)).^2).* clear_matrix;  
eye_raw_lowpass = eye_raw_lowpass(crop_x, :);
eye_visualize_lowpass = reshape(squeeze(eye_raw_lowpass(:,:)),[128 128]);  
[aa_lp, val_lp] = hist(eye_visualize_lowpass(:),linspace(0,max(...
                                    eye_visualize_lowpass(:)),1000));
    thresh_lp = val_lp(find(cumsum(aa_lp)/sum(aa_lp) > std_within,1,'first'));
eye_visualize_lowpass = uint16(eye_visualize_lowpass * 65536 / thresh_lp);

% For restored: linear interpolation
eye_raw_linear = sqrt(abs(squeeze(Data_img_linear1)).^2 + ...
           abs(squeeze(Data_img_linear2)).^2 + ...
           abs(squeeze(Data_img_linear3)).^2).* clear_matrix;  
eye_raw_linear = eye_raw_linear(crop_x, :);
eye_visualize_linear = reshape(squeeze(eye_raw_linear(:,:)),[128 128]);   
[aa_linear, val_linear] = hist(eye_visualize_linear(:),linspace(0,max(...
                                    eye_visualize_linear(:)),1000));
    thresh_linear = val_linear(find(cumsum(aa_linear)/sum(aa_linear) > std_within,1,'first'));
eye_visualize_linear = uint16(eye_visualize_linear * 65536 / thresh_linear);

% For restored: cubic interpolation
eye_raw_cubic = sqrt(abs(squeeze(Data_img_cubic1)).^2 + ...
           abs(squeeze(Data_img_cubic2)).^2 + ...
           abs(squeeze(Data_img_cubic3)).^2).* clear_matrix;  
eye_raw_cubic = eye_raw_cubic(crop_x, :);
eye_visualize_cubic = reshape(squeeze(eye_raw_cubic(:,:)),[128 128]);   
[aa_cubic, val_cubic] = hist(eye_visualize_cubic(:),linspace(0,max(...
                                    eye_visualize_cubic(:)),1000));
    thresh_cubic = val_cubic(find(cumsum(aa_cubic)/sum(aa_cubic) > std_within,1,'first'));
eye_visualize_cubic = uint16(eye_visualize_cubic * 65536 / thresh_cubic);

% For different channels (sub-images)    
eye_raw_ch1  = abs(squeeze(Data_img_good(:,:,1))).* clear_matrix;
eye_raw_ch2  = abs(squeeze(Data_img_good(:,:,2))).* clear_matrix;
eye_raw_ch3  = abs(squeeze(Data_img_good(:,:,3))).* clear_matrix;

eye_raw_ch1 = eye_raw_ch1(crop_x, :);
eye_raw_ch2 = eye_raw_ch2(crop_x, :);
eye_raw_ch3 = eye_raw_ch3(crop_x, :);

eye_visualize_ch1 = reshape(squeeze(eye_raw_ch1(:,:)),[128 128]); 
eye_visualize_ch2 = reshape(squeeze(eye_raw_ch2(:,:)),[128 128]); 
eye_visualize_ch3 = reshape(squeeze(eye_raw_ch3(:,:)),[128 128]); 

[aa_ch1, val_ch1] = hist(eye_visualize_ch1(:),linspace(0,max(...
                                    eye_visualize_ch1(:)),1000));
    thresh_ch1 = val_ch1(find(cumsum(aa_ch1)/sum(aa_ch1) > std_within,1,'first'));
    
[aa_ch2, val_ch2] = hist(eye_visualize_ch2(:),linspace(0,max(...
                                    eye_visualize_ch2(:)),1000));
    thresh_ch2 = val_ch2(find(cumsum(aa_ch2)/sum(aa_ch2) > std_within,1,'first'));
    
[aa_ch3, val_ch3] = hist(eye_visualize_ch3(:),linspace(0,max(...
                                    eye_visualize_ch3(:)),1000));
    thresh_ch3 = val_ch3(find(cumsum(aa_ch3)/sum(aa_ch3) > std_within,1,'first'));
    
eye_visualize_ch1 = uint16(eye_visualize_ch1 * 65536 / thresh_ch1); 
eye_visualize_ch2 = uint16(eye_visualize_ch2 * 65536 / thresh_ch2); 
eye_visualize_ch3 = uint16(eye_visualize_ch3 * 65536 / thresh_ch3);

% For restored: Wiener filter
eye_raw_wiener = sqrt(abs(squeeze(Data_img_wiener1)).^2 + ...
           abs(squeeze(Data_img_wiener2)).^2 + ...
           abs(squeeze(Data_img_wiener3)).^2).* clear_matrix;
eye_raw_wiener = eye_raw_wiener(crop_x, :);
eye_visualize_wiener = reshape(squeeze(eye_raw_wiener(:,:)),[128 128]);   
[aa_wiener, val_wiener] = hist(eye_visualize_wiener(:),linspace(0,max(...
                                    eye_visualize_wiener(:)),1000));
    thresh_wiener = val_wiener(find(cumsum(aa_wiener)/sum(aa_wiener) > std_within,1,'first'));
eye_visualize_wiener = uint16(eye_visualize_wiener * 65536 / thresh_wiener);

%% First fusion then restoration

eye_raw_bad_NEW = bad_data1 + bad_data2 + bad_data3 / 3;
eye_raw_good_NEW = good_data1 + good_data2 + good_data3 / 3;

difference = abs(eye_raw_good_NEW - eye_raw_bad_NEW);

for i = 1:size(difference,2)
    pow_col4(i) = sqrt(sum(difference(:,i).*difference(:,i)))/512;  
end

[TF4,P4] = islocalmax(pow_col4);
P4(TF4);

threshold4 = P4(6);

I4 = find(P4>=threshold4);

Inew4 = [I4(1:6),I4(8:end)];

%Inew4 = Inew1;

%imagesc(100*log(abs(eye_raw_good_NEW)));

% Linear Interpolation

linear_restored_new = eye_raw_bad_NEW;

for i = 1:length(Inew4)
    col1 = Inew4(i);
    linear_restored_new(:,col1) = (linear_restored_new(:,col1-1)+linear_restored_new(:,col1+1))/2;
end

Data_img_linear_new = ifftshift(ifft2(linear_restored_new),1);

eye_raw_linear_new = abs(squeeze(Data_img_linear_new)).* clear_matrix;
eye_raw_linear_new = eye_raw_linear_new(crop_x, :);
eye_visualize_linear_new = reshape(squeeze(eye_raw_linear_new(:,:)),[128 128]);   
[aa_ln, val_ln] = hist(eye_visualize_linear_new(:),linspace(0,max(...
                                    eye_visualize_linear_new(:)),1000));
    thresh_ln = val_ln(find(cumsum(aa_ln)/sum(aa_ln) > std_within,1,'first'));
eye_visualize_linear_new = uint16(eye_visualize_linear_new * 65536 / thresh_ln); 

% Cubic Interpolation

cubic_restored_new = eye_raw_bad_NEW;

for i = 1:length(Inew4)
    for j = 1:size(cubic_restored_new,1)
        coll3 = Inew4(i);
        seg3 = conv(cubic_filt,cubic_restored_new(j,coll3-5:coll3+5),'same');
        cubic_restored_new(j,coll3) = seg3(4);
    end
end

Data_img_cubic_new = ifftshift(ifft2(cubic_restored_new),1);

eye_raw_cubic_new = abs(squeeze(Data_img_cubic_new)).* clear_matrix;
eye_raw_cubic_new = eye_raw_cubic_new(crop_x, :);
eye_visualize_cubic_new = reshape(squeeze(eye_raw_cubic_new(:,:)),[128 128]);   
[aa_cn, val_cn] = hist(eye_visualize_cubic_new(:),linspace(0,max(...
                                    eye_visualize_cubic_new(:)),1000));
    thresh_cn = val_cn(find(cumsum(aa_cn)/sum(aa_cn) > std_within,1,'first'));
eye_visualize_cubic_new = uint16(eye_raw_linear_new * 65536 / thresh_cn); 

% Gaussian Filter

gaussian_restored_new = eye_raw_bad_NEW;

for i = 1:length(Inew4)
    colll1 = Inew4(i);
    segg1 = conv2(gaussian_restored_new(:,colll1-1:colll1+1),gaus_filter,'same');
    gaussian_restored_new(:,colll1) = segg1(:,2);
end

Data_img_gauss_new = ifftshift(ifft2(gaussian_restored_new),1);

eye_raw_gaussian_new  = abs(squeeze(Data_img_gauss_new)).* clear_matrix;  
eye_raw_gaussian_new = eye_raw_gaussian_new(crop_x, :);
eye_visualize_gaussian_new = reshape(squeeze(eye_raw_gaussian_new(:,:)),[128 128]);  
[aa_gauss, val_gauss] = hist(eye_visualize_gaussian_new(:),linspace(0,max(...
                                    eye_visualize_gaussian_new(:)),1000));
    thresh_gauss = val_gauss(find(cumsum(aa_gauss)/sum(aa_gauss) > std_within,1,'first'));
eye_visualize_gaussian_new = uint16(eye_visualize_gaussian_new * 65536 / thresh_gauss);

% Low-Pass Filter

lowpass_restored_new = eye_raw_bad_NEW;

for i = 1:length(Inew4)
    colll1 = Inew4(i);
    segg1 = conv2(lowpass_restored_new(:,colll1-1:colll1+1),lowpass_filter,'same');
    lowpass_restored_new(:,colll1) = segg1(:,2);
end

Data_img_lowpass_new = ifftshift(ifft2(lowpass_restored_new),1);

eye_raw_lowpass_new  = abs(squeeze(Data_img_lowpass_new)).* clear_matrix;  
eye_raw_lowpass_new = eye_raw_lowpass_new(crop_x, :);
eye_visualize_lowpass_new = reshape(squeeze(eye_raw_lowpass_new(:,:)),[128 128]);  
[aa_lp, val_lp] = hist(eye_visualize_lowpass_new(:),linspace(0,max(...
                                    eye_visualize_lowpass_new(:)),1000));
    thresh_lp = val_lp(find(cumsum(aa_lp)/sum(aa_lp) > std_within,1,'first'));
eye_visualize_lowpass_new = uint16(eye_visualize_lowpass_new * 65536 / thresh_lp);

% Wiener Filter

wiener_restored_new = eye_raw_bad_NEW;

for i = 1:length(Inew4)
    colll1 = Inew4(i);
    segg1 = conv2(eye_raw_bad_NEW(:,colll1-nn:colll1+nn),filter1,'same');
    wiener_restored_new(:,colll1) = segg1(:,2);
end

Data_img_wiener_new = ifftshift(ifft2(wiener_restored_new),1);

eye_raw_wiener_new = abs(squeeze(Data_img_wiener_new)).* clear_matrix;
eye_raw_wiener_new = eye_raw_wiener_new(crop_x, :);
eye_visualize_wiener_new = reshape(squeeze(eye_raw_wiener_new(:,:)),[128 128]);   
[aa_wiener, val_wiener] = hist(eye_visualize_wiener_new(:),linspace(0,max(...
                                    eye_visualize_wiener_new(:)),1000));
    thresh_wiener = val_wiener(find(cumsum(aa_wiener)/sum(aa_wiener) > std_within,1,'first'));
eye_visualize_wiener_new = uint16(eye_visualize_wiener_new * 65536 / thresh_wiener);

%% Variance calculation for three sub-channels
region1 = eye_visualize_ch1(110:115,110:115);
region2 = eye_visualize_ch2(110:115,110:115);
region3 = eye_visualize_ch3(110:115,110:115);

var1 = VarianceFinder(region1);
var2 = VarianceFinder(region2);
var3 = VarianceFinder(region3);

mean = (var1+var2+var3)/3;

%% SNR comparisons

uncropped_good = good_data1;
uncropped_bad = bad_data1;
uncropped_gauss = gauss_restored1;
uncropped_lowpass = lowpass_restored1;
uncropped_linear = linear_restored1;
uncropped_cubic = cubic_restored1;
uncropped_wiener = wiener_restored1;

uncropped_gaussian_new = gaussian_restored_new;
uncropped_lowpass_new = lowpass_restored_new;
uncropped_linear_new = linear_restored_new;
uncropped_cubic_new = cubic_restored_new;
uncropped_wiener_new = wiener_restored_new;

%-

uncropped_good_spatial = ifftshift(ifft2(uncropped_good),1);
uncropped_bad_spatial = ifftshift(ifft2(uncropped_bad),1);
uncropped_gauss_spatial = ifftshift(ifft2(uncropped_gauss),1);
uncropped_lowpass_spatial = ifftshift(ifft2(uncropped_lowpass),1);
uncropped_linear_spatial = ifftshift(ifft2(uncropped_linear),1);
uncropped_cubic_spatial = ifftshift(ifft2(uncropped_cubic),1);
uncropped_wiener_spatial = ifftshift(ifft2(uncropped_wiener),1);

uncropped_gaussain_spatial_new = ifftshift(ifft2(uncropped_gaussian_new),1);
uncropped_lowpass_spatial_new = ifftshift(ifft2(uncropped_lowpass_new),1);
uncropped_linear_spatial_new = ifftshift(ifft2(uncropped_linear_new),1);
uncropped_cubic_spatial_new = ifftshift(ifft2(uncropped_cubic_new),1);
uncropped_wiener_spatial_new = ifftshift(ifft2(uncropped_wiener_new),1);

%-

noisy_region = uncropped_bad_spatial(10:70,:);
var_noise = VarianceFinder(noisy_region);


noisy_bad_fused = ifftshift(ifft2(eye_raw_bad_NEW),1);
noisy_region_new = noisy_bad_fused(10:70,:);

var_noise_new = VarianceFinder(noisy_region_new);

var_bad = VarianceFinder(Data_img_bad(:,:,1));
var_gauss = VarianceFinder(Data_img_gauss1);
var_lowpass = VarianceFinder(Data_img_lowpass1);
var_linear = VarianceFinder(Data_img_linear1);
var_cubic = VarianceFinder(Data_img_cubic1);
var_wiener = VarianceFinder(Data_img_wiener1);

var_gauss_new = VarianceFinder(Data_img_gauss_new);
var_lowpass_new = VarianceFinder(Data_img_lowpass_new);
var_linear_new = VarianceFinder(Data_img_linear_new);
var_cubic_new = VarianceFinder(Data_img_cubic_new);
var_wiener_new = VarianceFinder(Data_img_wiener_new);


SNR_bad = SNRfinder(var_bad,var_noise);
SNR_gauss = SNRfinder(var_gauss,var_noise);
SNR_lowpass = SNRfinder(var_lowpass,var_noise);
SNR_linear = SNRfinder(var_linear,var_noise);
SNR_cubic = SNRfinder(var_cubic,var_noise);
SNR_wiener = SNRfinder(var_wiener,var_noise);

SNR_gauss_new = SNRfinder(var_gauss_new,var_noise_new);
SNR_lowpass_new = SNRfinder(var_lowpass_new,var_noise_new);
SNR_linear_new = SNRfinder(var_linear_new,var_noise_new);
SNR_cubic_new = SNRfinder(var_cubic_new,var_noise_new);
SNR_wiener_new = SNRfinder(var_wiener_new,var_noise_new);


fprintf('SNR of Original Bad Image = %2.5f dB\n',SNR_bad);

fprintf('Restoration then Fusion\n');
fprintf('SNR of Gaussian Filter Image = %2.5f dB\n',SNR_gauss);
fprintf('SNR of Low-Pass Filter Image = %2.5f dB\n',SNR_lowpass);
fprintf('SNR of Linear Interpolation Image = %2.5f dB\n',SNR_linear);
fprintf('SNR of Cubic Interpolation Image = %2.5f dB\n',SNR_cubic);
fprintf('SNR of Wiener Filter Image = %2.5f dB\n',SNR_wiener);

fprintf('Fusion then Restoration\n');
fprintf('SNR of Gaussian Filter Image = %2.5f dB\n',SNR_gauss_new);
fprintf('SNR of Low-Pass Filter Image = %2.5f dB\n',SNR_lowpass_new);
fprintf('SNR of Linear Interpolation Image = %2.5f dB\n',SNR_linear_new);
fprintf('SNR of Cubic Interpolation Image = %2.5f dB\n',SNR_cubic_new);
fprintf('SNR of Wiener Filter Image = %2.5f dB\n',SNR_wiener_new);


%% Optimal Fusion Problem

% For original: good data with optimally fusing
eye_raw_good_optimal = sqrt((mean/var1)*abs(squeeze(Data_img_good(:,:,1))).^2 + ...
           abs((mean/var2)*squeeze(Data_img_good(:,:,2))).^2 + ...
           abs((mean/var3)*squeeze(Data_img_good(:,:,3))).^2).* clear_matrix;  
eye_raw_good_optimal = eye_raw_good_optimal(crop_x, :);
eye_visualize_good_optimal = reshape(squeeze(eye_raw_good_optimal(:,:)),[128 128]);   
[aa_good_optimal, val_good_optimal] = hist(eye_visualize_good_optimal(:),linspace(0,max(...
                                    eye_visualize_good_optimal(:)),1000));
    thresh_good_optimal = val_good_optimal(find(cumsum(aa_good_optimal)/sum(aa_good_optimal) > std_within,1,'first'));
eye_visualize_good_optimal = uint16(eye_visualize_good_optimal * 65536 / thresh_good_optimal); 

%% Plotting

figure(1)
imagesc(eye_visualize_good);
title("Fused Image, Slice " + slice_number + "(Good Data)")
axis image, 
colormap gray;
axis off

figure(2); 
sgtitle("Channel Images, Slice " + slice_number + "(Good Data)")
subplot(1,3,1)
imagesc(eye_visualize_ch1);
axis image, 
colormap gray;
axis off
title("Channel 1")
subplot(1,3,2)
imagesc(eye_visualize_ch2);
axis image, 
colormap gray;
axis off
title("Channel 2")
subplot(1,3,3)
imagesc(eye_visualize_ch3);
title("Channel 3")
axis image, 
colormap gray;
axis off

figure(3); 
sgtitle("Spatial Frequency Observations, Slice " + slice_number + "(Good Data)")
subplot(1,3,1)
imagesc(100*log(abs(good_data1)));
title("Channel 1")
xlabel('Horizontal frequency bins')
ylabel('Vertical frequency bins');
subplot(1,3,2)
imagesc(100*log(abs(good_data2)));
title("Channel 2")
xlabel('Horizontal frequency bins')
ylabel('Vertical frequency bins');
subplot(1,3,3)
imagesc(100*log(abs(good_data3)));
title("Channel 3")
xlabel('Horizontal frequency bins')
ylabel('Vertical frequency bins');

figure(4); 
sgtitle("Spatial Frequency Observations, Slice " + slice_number + "(Bad Data)")
subplot(1,3,1)
imagesc(100*log(abs(bad_data1)));
title("Channel 1")
xlabel('Horizontal frequency bins')
ylabel('Vertical frequency bins');
subplot(1,3,2)
imagesc(100*log(abs(bad_data2)));
title("Channel 2")
xlabel('Horizontal frequency bins')
ylabel('Vertical frequency bins');
subplot(1,3,3)
imagesc(100*log(abs(bad_data3)));
title("Channel 3")
xlabel('Horizontal frequency bins')
ylabel('Vertical frequency bins');

figure(5);
sgtitle("Outlier Detection and Rejection, Slice " + slice_number + "Channel 1 (Bad Data)")
subplot(1,5,1)
imagesc(100*log(abs(bad_data1)));
title("With Outliers")
xlabel('Horizontal frequency bins')
ylabel('Vertical frequency bins');
subplot(1,5,2)
imagesc(100*log(abs(outlier_data1)));
title("Outliers Detected")
xlabel('Horizontal frequency bins')
ylabel('Vertical frequency bins');
subplot(1,5,3)
imagesc(100*log(abs(linear_restored1)));
title("Outliers Rejected with Linear Interpolation")
xlabel('Horizontal frequency bins')
ylabel('Vertical frequency bins');
subplot(1,5,4)
imagesc(100*log(abs(cubic_restored1)));
title("Outliers Rejected with Cubic Interpolation")
xlabel('Horizontal frequency bins')
ylabel('Vertical frequency bins');

figure(6)
sgtitle("Fused Images, Slice " + slice_number + "(Bad Data)")
subplot(1,3,1)
imagesc(eye_visualize_bad)
axis image, 
colormap gray;
axis off
title('Original')
subplot(1,3,2)
imagesc(eye_visualize_linear)
axis image, 
colormap gray;
axis off
title('Linear Interpolation')
subplot(1,3,3)
imagesc(eye_visualize_cubic)
axis image, 
colormap gray;
axis off
title('Cubic Interpolation')

figure(7)
sgtitle("Fused Images, Slice " + slice_number + "(Bad Data)")
subplot(1,2,1)
imagesc(eye_visualize_bad)
axis image, 
colormap gray;
axis off
title('Original')
subplot(1,2,2)
imagesc(eye_visualize_gaussian)
axis image, 
colormap gray;
axis off
title('With a Gaussian Filter')

figure(8)
sgtitle("Fused Images, Slice " + slice_number + "(Bad Data)")
subplot(1,2,1)
imagesc(eye_visualize_bad)
axis image, 
colormap gray;
axis off
title('Original')
subplot(1,2,2)
imagesc(eye_visualize_lowpass)
axis image, 
colormap gray;
axis off
title('With a Low-Pass Filter')

figure(9)
sgtitle("Restoration and Fusion Order, Slice " + slice_number + "(Bad Data)")

subplot(2,5,1)
imagesc(eye_visualize_lowpass)
axis image, 
colormap gray;
axis off
title('Low-Pass (Restore, Fuse)')

subplot(2,5,2)
imagesc(eye_visualize_gaussian)
axis image, 
colormap gray;
axis off
title('Gauss (Restore, Fuse)')

subplot(2,5,3)
imagesc(eye_visualize_linear)
axis image, 
colormap gray;
axis off
title('Linear (Restore, Fuse)')

subplot(2,5,4)
imagesc(eye_visualize_cubic)
axis image, 
colormap gray;
axis off
title('Cubic (Restore, Fuse)')

subplot(2,5,5)
imagesc(eye_visualize_wiener)
axis image, 
colormap gray;
axis off
title('Wiener (Restore, Fuse)')

subplot(2,5,6)
imagesc(eye_visualize_lowpass_new)
axis image, 
colormap gray;
axis off
title('Low-Pass (Fuse, Restore)')

subplot(2,5,7)
imagesc(eye_visualize_gaussian_new)
axis image, 
colormap gray;
axis off
title('Gauss (Fuse, Restore)')

subplot(2,5,8)
imagesc(eye_visualize_linear_new)
axis image, 
colormap gray;
axis off
title('Linear (Fuse, Restore)')

subplot(2,5,9)
imagesc(eye_visualize_cubic_new)
axis image, 
colormap gray;
axis off
title('Cubic (Fuse, Restore)')

subplot(2,5,10)
imagesc(eye_visualize_wiener_new)
axis image, 
colormap gray;
axis off
title('Wiener (Fuse, Restore)')

figure(10);
sgtitle("Outlier Detection, Slice " + slice_number + ", Channel 1 (Bad Data)")
subplot(1,2,1)
imagesc(100*log(abs(bad_data1)));
title("Outliers")
xlabel('Horizontal frequency bins')
ylabel('Vertical frequency bins');
subplot(1,2,2)
imagesc(100*log(abs(outlier_data1)));
title("Outliers Detected")
xlabel('Horizontal frequency bins')
ylabel('Vertical frequency bins');

figure(11)
sgtitle("Fused Images, Slice " + slice_number + "(Bad Data)")
subplot(1,2,1)
imagesc(eye_visualize_bad)
axis image, 
colormap gray;
axis off
title('Original')
subplot(1,2,2)
imagesc(eye_visualize_wiener)
axis image, 
colormap gray;
axis off
title('With a Wiener Filter')

figure(12)
subplot(1,2,1)
sgtitle("Fused Images, Slice " + slice_number + "(Good Data)")
imagesc(eye_visualize_good);
title("Original")
axis image, 
colormap gray;
axis off
subplot(1,2,2)
sgtitle("Fused Images, Slice " + slice_number + "(Good Data)")
imagesc(eye_visualize_good_optimal);
title("Optimal Fusing")
axis image, 
colormap gray;
axis off

%% Functions

function varian_square = VarianceFinder(region)

    %Type conversion
    region = double(region);

    %Getting the size of the region
    N1 = size(region,1);
    N2 = size(region,2);

    %Finding the average of the region
    avg = sum(region,'all')/(N1*N2);

    summation = 0;

    %Iterating through the region for finding the variance V
    for i = 1:N1
        for j = 1:N2
            summation = summation+(region(i,j)-avg)^2; 
        end
    end

    %Finding the variance V of the flat region
    varian_square = summation/(N1*N2);

end

function SNR = SNRfinder(var_signal,var_noise)
    SNR = 10*log10(var_signal/var_noise);
end
