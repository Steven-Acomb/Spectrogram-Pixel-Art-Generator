clear all %#ok<CLALL> 
clc

% Author: Stephen Acomb
% License: GPLv3
% Maintainer: Stephen Acomb
% Email: acomb.stephen@gmail.com

% Define colormap of 255 hues plus black
% uniform_h = linspace(0,1,255);
unshifted_hues = [linspace(0,1,255) linspace(0,1,255)];
shift = 0;
uniform_h = NaN(1,255);
for k = 1:255
    uniform_h(k) = unshifted_hues(k+shift);
end
uniform_hsv = NaN(255,1,3);
for k = 1:255
    uniform_hsv(k,1,:) = [uniform_h(k) 1 1];
end
uniform_rgb = hsv2rgb(uniform_hsv);
% imshow(uniform_rgb)
uniform = NaN(256,3);
for k = 2:256
    uniform(k,:) = uniform_rgb(k-1,1,:);
end
uniform(1,:) = [0, 0, 0];

% Input Variables:
N = 64; % Image resolution is NxN pixels.
xRng = [0 10];
yRng = [0 10];
td = 1;
fs = 4e4;
fMax = 0.25*fs;
fMin = 0;
dBpk = -1;
% specColormap = "turbo";
specColormap = uniform;
specFloor_dB = -63;
% specFloor_dB = -45;
% specFloor_dB = -60;
% specFloor_dB = -90;
% specFloor_dB = -Inf;
specOneAtATime = false;
logf = false;
[filenameOut,filepathOut] = uiputfile('*.wav');
audioFilename = fullfile(filepathOut,filenameOut);
% audioFilename = 'C:\Users\acombsr\OneDrive - Rose-Hulman Institute of Technology\Desktop\sound.wav';

% Define Constants:
del = N*2^-10;

% imgin = procureInputImage(N);
% imshow(imgin)

% Prompt the user for a PNG & convert it to an NxN array of colors that
% are either black (blank on spectrogram) or vary only by their hue
% in order to reduce the number of dimensions.
% function output_image = procureInputImage(N)
    % Get input image
%     file = 'C:\Users\acombsr\OneDrive - Rose-Hulman Institute of Technology\Desktop\Work\481\Mini-Projects\MP5 - Additive Synthesis\test images\rainbow.jpg';
    [filename,filepath] = uigetfile({'*.png';'*.jpg'});
    file = fullfile(filepath,filename);
    input_image = imread(file);
    
    % Ensure image is scaled correctly, then convert to HSV color
    rgb_image = imresize(input_image,[N N]);
    hsv_image = rgb2hsv(rgb_image);
    
%     adj_image = hsv_image;
    hue_map = NaN(N,N);
    for i = 1:size(hsv_image,1)
        for j = 1:size(hsv_image,2)
    %         if (hsv_image(i,j,3) < 0.5)  (hsv_image(i,j,2) < 0.5)
            if hsv_image(i,j,3) >= 0.5
                hue_map(i,j) = hsv_image(i,j);
            end
%             if hsv_image(i,j,3) < 0.5
%                 adj_image(i,j,3) = 0;
%             else
%                 adj_image(i,j,3) = 1;
%             end
%             adj_image(i,j,2) = 1;
        end
    end
%     output_image = hsv2rgb(adj_image);
% end

% % Reverse conversion here for debugging purposes
adj_image = NaN(N,N,3);
for i = 1:size(hsv_image,1)
    for j = 1:size(hsv_image,2)
        if isnan(hue_map(i,j))
            adj_image(i,j,1) = 0;
            adj_image(i,j,3) = 0;
        else
            adj_image(i,j,1) = hue_map(i,j);
            adj_image(i,j,3) = 1;
        end
        adj_image(i,j,2) = 1;
    end
end
output_image = hsv2rgb(adj_image);
converted_file = 'C:\Users\acombsr\OneDrive - Rose-Hulman Institute of Technology\Desktop\Work\481\Mini-Projects\MP5 - Additive Synthesis\test images\rainbow_converted.jpg';
imwrite(output_image, converted_file)
% imshow(output_image)

% Create Partials
% partials = repmat(TimbrePartial,1,N+2);
partials = repmat(TimbrePartial,1,N);
for p = 1:N
    n = hue_map(p,:);
    partials(p) = pixelRowPartial(n,1+N-p,N,xRng);
end

% partials(N+1).fy = [0.5, 0.5];
% partials(N+1).fx = xRng;
% partials(N+1).fim = 'linear';
% partials(N+1).ay = [1, 1];
% partials(N+1).ax = xRng;
% partials(N+1).aim = 'linear';
% 
% partials(N+2).fy = [17.5, 17.5];
% partials(N+2).fx = xRng;
% partials(N+2).fim = 'linear';
% partials(N+2).ay = [1, 1];
% partials(N+2).ax = xRng;
% partials(N+2).aim = 'linear';

% p = 2;
% partials = repmat(pixelRowPartial(hue_map(p,:),p,N,xRng),1,1);

% function partial = pixelRowPartial(n,k,N,xRng)

% p = 1;
% n = hue_map(p,:);
% % n = linspace(0,1,N);
% k = p;
% 
% % function partial = pixelRowPartial(n,k,N,xRng)
%     minval = 10^(-1*((digits/2)-2));
%     xL = xRng(1);
%     xH = xRng(2);
%     del = N * 2^-10;
%     fy = [1+del*k, 1+del*k];
%     fx = [xL, xH];
%     fim = 'linear';
%     anchors = linspace(xL,xH,N+3);
%     ax = zeros(1,3*(N+1));
%     ay = (-Inf)*ones(1,3*(N+1));
%     for i = 1:3*(N+1)
%         if mod((i),3) == 2
%             ax(i) = anchors(1+(i+1)/3);
%         end
%         switch mod((i),3)
%             case 0
%                 ax(i) = anchors(1+floor((i+1)/3)) + minval;
%             case 1
%                 ax(i) = anchors(1+floor((i+1)/3)) - minval;
%             case 2
%                 ax(i) = anchors(1+(i+1)/3);
%         end
%     end
%     for j = 1:N
%         switch n(j)
%             case 0
%                 ay(3*j) = -99;
%                 ay(3*j+1) = -99;
%             otherwise
%                 ay(3*j) = 10*n(j);
%                 ay(3*j+1) = 10*n(j);
%         end    
%     end
%     ax(1) = xL;
%     ax(3*(N+1)) = xH;
%     aim = 'linear';
%     partial = TimbrePartial(fy,fx,fim,ay,ax,aim);
% % end
% 
% partials = partial;

% partial_image = NaN(1,N,3);
% for j = 1:N
%     if isnan(hue_map(i,j))
%         adj_image(i,j,1) = 0;
%         adj_image(i,j,3) = 0;
%     else
%         adj_image(i,j,1) = hue_map(i,j);
%         adj_image(i,j,3) = 1;
%     end
%     adj_image(i,j,2) = 1;
% end

% Create Partials
patch = AdditiveSynthPatch(partials,xRng,yRng,fs,td,logf,...
    fMax,fMin,dBpk);
% patch.show(specFloor_dB,specColormap);
% patch.show_freq_traj();
% patch.show_ampl_traj();
patch.show_spec(specFloor_dB,specColormap);
patch.listen(specOneAtATime);
patch.write(audioFilename);


function partial = pixelRowPartial(n,k,N,xRng)
    minval = 10^(-1*((digits/2)-2));
    xL = xRng(1);
    xH = xRng(2);
    del = 16/N;
    fy = [1+del*k, 1+del*k];
    fx = [xL, xH];
    fim = 'linear';
    anchors = linspace(xL,xH,N+3);
    ax = zeros(1,3*(N+1));
    ay = (-Inf)*ones(1,3*(N+1));
    for i = 1:3*(N+1)
        if mod((i),3) == 2
            ax(i) = anchors(1+(i+1)/3);
        end
        switch mod((i),3)
            case 0
                ax(i) = anchors(1+floor((i+1)/3)) + minval;
            case 1
                ax(i) = anchors(1+floor((i+1)/3)) - minval;
            case 2
                ax(i) = anchors(1+(i+1)/3);
        end
    end
    for j = 1:N
        switch n(j)
            case 0
                ay(3*j) = -99;
                ay(3*j+1) = -99;
            otherwise
                ay(3*j) = 10*n(j);
                ay(3*j+1) = 10*n(j);
        end    
    end
    ax(1) = xL;
    ax(3*(N+1)) = xH;
    aim = 'linear';
    partial = TimbrePartial(fy,fx,fim,ay,ax,aim);
end