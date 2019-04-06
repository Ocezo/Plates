% Plate Recognition - 03/10/2019
% Jean-Marc Berthommé
%
% - 11/27/2014:
%   . 1st version: plate crop from feature extraction
% - 12/02/2014:
%   . Color clustering more or less initialized
% - 12/04/2014:
%   . Color labelling improved (1st pass removal + special k-means)
%   . Harris corners added
% - 05/20/2017:
%   . reworking on the display
% - 03/10/2019:
%   . Corners detection update (see git commit message)

function plate_recognition
dbg = 1;      % debug flag
i = 8;        % image index - 1-35
nam = 'car';  % image name
dir = 'Img';  % image directory
ext = 'jpg';  % image extension

% read and display the source image
file = sprintf('%s%d.%s',nam,i,ext);
path = sprintf('%s/%s',dir,file);
I = imread(path);
disp_src(I,file);

% graphical input
if     i== 4, x = [600; 1600]; y = [452; 784];  % <~  car4.jpg
elseif i== 5, x = [776; 1162]; y = [534; 684];  % <~  car5.jpg
elseif i== 6, x = [752; 1096]; y = [484; 686];  % <~  car6.jpg
elseif i== 7, x = [663; 1199]; y = [530; 684];  % <~  car7.jpg
elseif i== 8, x = [713; 1486]; y = [348; 644];  % <~  car8.jpg
elseif i== 9, x = [831; 1138]; y = [371; 634];  % <~  car9.jpg
elseif i==10, x = [900; 1295]; y = [744; 880];  % <~ car10.jpg
elseif i==11, x = [777; 1418]; y = [775; 1007]; % <~ car11.jpg
elseif i==12, x = [636; 1368]; y = [598;  838]; % <~ car12.jpg
else         [x,y] = ginput(2);
end

% region of interest (ROI)
y = round(y); x = round(x);
h = y(2)-y(1)+1; w = x(2)-x(1)+1;
disp_rect([x(1),y(1)],h,w,'y',2);

% crop the image and display it
Ic = I(y(1):y(2),x(1):x(2),:);
disp_crop(Ic);

% expand the RGB data
X = reshape(double(Ic), h*w, 3);

% main plate colors
mpcol  = mean(X,1);     % mean plate color
white  = [250 240 235];
black  = [ 95  85  80];
dblue  = [ 75 135 185]; % dark  blue
lblue  = [130 255 255]; % light blue
red    = [235   5  10];
yellow = [210 160  10];
pcols  = [white; black; dblue; lblue; red; yellow]; % desired plate colors
disp_col(mpcol,pcols);

% get all the colors distances
d = euc_dist(X,pcols);

% display the threshed maps of the desired plate colors
th = 60;
disp_lab(d,th,h,w);

% eliminate some plate colors
pct = sum(d < th, 1) / (h*w);  % percentages list
p = (pct > 0.01)';             % "present" colors
pcols = pcols(p,:);

% kmeans with a color initialization
opt = 'col';                   % color option: art / col
k = 1 + sum(p);                % total number of clusters
C = kmeans(X,[mpcol;pcols],1); Imap = map(opt,k,X,C,h,w);

% image corners (from Harris & Stephen or Noble)
Ig = rgb2gray(Ic); % Grayscale conversion
thresh = 20;       % Harris threshold \in [10, 50]
sigma = 1;         % dilution parameter
radius = 2;        % threshold radius
% kh = 0.04;       % Harris's "k" (not to confuse with k-means's "k")
% [~,r,c] = harris(Ig, sigma, kh, 'thresh', thresh, 'radius', radius);
[~,r,c] = noble(Ig, sigma, 'thresh', thresh, 'radius', radius);

Hc=[r,c]; % Harris corners
Ih = encrust_harris(Imap,Hc); disp_har(dbg,Ih,Hc);

% Absolute gradient images
[Ix, Iy] = derivative5(Ig,'x','y');
Iax = abs(Ix); Iay = abs(Iy);
disp_abs_grad(dbg, Iax, Iay);

% Gradient norm and angle images
Nxy = (Ix.^2 + Iy.^2).^0.5; % gradient norm
Axy = atan2(-Iy,-Ix);       % gradient angles on [-pi, pi] with convention:
                            % .--> x  
                            % |
                            % v y

% This code aims to detect fast the horizontal pararel lines of the plates
% in order to ameliorate the 2D homography
% *********   *********
[NXY, ~] = disp_norm_angle(dbg, Nxy, Axy);

set(7,'Position',[1200 330 720 645]);

fprintf('[min(Nxy), max(Nxy)]= [%0.2f, %0.2f]\n',min(Nxy(:)),max(Nxy(:)));
fprintf('[min(NXY), max(NXY)]= [%0.2f, %0.2f]\n',min(NXY(:)),max(NXY(:)));

H = histc(NXY(:), 0:255); % gray-level histogram of the Norm of image
figure; bar(H);

P = H/(h*w);    % sorted pdf - from 0 to 255
CP = cumsum(P); % sorted cdf - from 0 to 255

% sum(P)        % check that \Sigma p_k = 1 !
% [P, CP]       % display the pdf & the cdf ;D

Test = CP > 0.95;    % get 95% of the darkest pixels (~ less +/- 3 sigma)

UINT8 = (0:255)';    % vector of gray-level 
[P, CP, Test, UINT8]

% find teh the threshold
TH = min(UINT8(Test));

NXY_TH = NXY > TH; % build the threshed image - 84
figure; image(repmat(uint8(255*NXY_TH), [1 1 3])); axis image;
% **************************

% find the main horizontal //
% ********** TODO **********
Dxy = Nxy > (mean(Nxy(:)) + std(Nxy(:)));
figure; image(repmat(Dxy, [1 1 3])); axis image;

Dxy = Nxy > (mean(Nxy(:)) + 2*std(Nxy(:)));
figure; image(repmat(Dxy, [1 1 3])); axis image;

Dxy = Nxy > (mean(Nxy(:)) + 3*std(Nxy(:)));
figure; image(repmat(Dxy, [1 1 3])); axis image;

A = Axy(Dxy); W = Iay(Dxy);
% figure; hist(A, linspace(0,pi,1000));
bins = 10000; amin = min(A); amax = max(A);
[histw, histv] = histwv(A, W, amin, amax, bins);

disp_hist(dbg, histw, histv);

delta = (amax-amin)/(bins-1); [~,idmax] = max(histw);
theta_max = (idmax-1)*delta + amin;
fprintf('L''angle max est de %0.2f°.\n\n', theta_max * 180/pi);
% **************************

% Canny edges
Ie = edge(Ig,'canny');
disp_can(dbg, Ie);

% Igmap = rgb2gray(Imap);
% Iemap = edge(Igmap,'canny');
% figure; image(repmat(Iemap,[1 1 3])); axis image; axis off;

% descriptors

% 2d Homography based on good descriptors

% visual parsing

% call tesseract or gocr

fprintf('Press any key to continue...\n');
pause; clear all; close all;

function disp_src(I, file)
% display the source image
tit = sprintf('Raw Image - %s', file);
f1 = figure(1); set(gcf,'Color',[0.2,0.2,0.2]);
set(f1,'Position',[1922 552 635 445]);
image(I); title(tit,'Color','w'); axis off; axis image;

function disp_rect(X, height, width, color, linewidth)
% display a rectangle that enlights the ROI
Px = [X(1); X(1)+width-1; X(1)+width-1; X(1); X(1)];
Py = [X(2); X(2); X(2)+height-1; X(2)+height-1; X(2)];
hold on; plot(Px,Py,'Color',color, 'LineWidth', linewidth); hold off;

function disp_crop(Ic)
% display the cropped image
f2 = figure(2); set(gcf,'Color',[0.2,0.2,0.2]);
set(f2,'Position',[2570 556 635 445]);
image(Ic); title('Crop','Color','w'); axis off; axis image;

function disp_col(mpcol, pcols)
% display the desired plate colors
f3 = figure(3); set(gcf,'Color',[0.2,0.2,0.2]);
set(f3,'Position',[1926 6 635 445]);
subplot(2,4,1); image(reshape(uint8(mpcol), [1 1 3]));
title('Mean Plate Color','Color','w'); axis off; axis image;
subplot(2,4,3); image(reshape(uint8(pcols(1,:)),[1 1 3]));
title('White','Color','w'); axis off; axis image;
subplot(2,4,4); image(reshape(uint8(pcols(2,:)),[1 1 3]));
title('Black','Color','w'); axis off; axis image;
subplot(2,4,5); image(reshape(uint8(pcols(3,:)), [1 1 3]));
title('Dark Blue','Color','w'); axis off; axis image;
subplot(2,4,6); image(reshape(uint8(pcols(4,:)), [1 1 3]));
title('Light Blue','Color','w'); axis off; axis image;
subplot(2,4,7); image(reshape(uint8(pcols(5,:)), [1 1 3]));
title('Red','Color','w'); axis off; axis image;
subplot(2,4,8); image(reshape(uint8(pcols(6,:)), [1 1 3]));
title('Yellow','Color','w'); axis off; axis image;

function d = euc_dist(X,Y)
% calculate all the Euclidean distances between 2 populations in d-space
nx = size(X,1); ny = size(Y,1);
d = sqrt(sum(X.^2,2)*ones(1,ny) + ones(nx,1)*sum(Y.^2,2)' - 2*(X*Y'));

function disp_lab(d, th, h, w)
% display labels
% [0        , (255^2*3)^0.5]
% [min(d(:)),     max(d(:))]
% figure; hist(d(:), linspace(0,441.68,50));
dt = d < th;                % threshed color distances
pct = sum(dt,1) / (h*w);    % area pct
COL = [0 0 0; 255 255 255]; % black & white

Cw = dt(:,1); Iw = reshape(uint8(COL(Cw+1,:)),h,w,3); % white
Cb = dt(:,2); Ib = reshape(uint8(COL(Cb+1,:)),h,w,3); % black
Cd = dt(:,3); Id = reshape(uint8(COL(Cd+1,:)),h,w,3); % dark blue
Cl = dt(:,4); Il = reshape(uint8(COL(Cl+1,:)),h,w,3); % light blue
Cr = dt(:,5); Ir = reshape(uint8(COL(Cr+1,:)),h,w,3); % red
Cy = dt(:,6); Iy = reshape(uint8(COL(Cy+1,:)),h,w,3); % yellow

f4 = figure(4); set(gcf,'Color',[0.2,0.2,0.2]);
set(f4,'Position',[2571 10 635 445]);

titw = sprintf('White - th = %d -> %0.1f %%',th, pct(1)*100);
titb = sprintf('Black - th = %d -> %0.1f %%',th, pct(2)*100);
titd = sprintf('Dark Blue - th = %d -> %0.1f %%',th, pct(3)*100);
titl = sprintf('Light Blue - th = %d -> %0.1f %%',th, pct(4)*100);
titr = sprintf('Red - th = %d -> %0.1f %%',th, pct(5)*100);
tity = sprintf('Yellow - th = %d -> %0.1f %%',th, pct(6)*100);

subplot(2,3,1);image(Iw); title(titw,'Color','w'); axis off; axis image;
subplot(2,3,2);image(Ib); title(titb,'Color','w'); axis off; axis image;
subplot(2,3,3);image(Id); title(titd,'Color','w'); axis off; axis image;
subplot(2,3,4);image(Il); title(titl,'Color','w'); axis off; axis image;
subplot(2,3,5);image(Ir); title(titr,'Color','w'); axis off; axis image;
subplot(2,3,6);image(Iy); title(tity,'Color','w'); axis off; axis image;

function Imap = map(opt, k, X, C, h, w)
% display the image mapping
if strcmp(opt, 'art')
    COL = 255*hsv(k);  % HSV colors
elseif strcmp(opt, 'col')
    COL = zeros(k, 3); % RGB colors
    for i=1:k, COL(i,:) = uint8(mean(X(C==i,:),1)); end;
else
    error('Unknown "%s" color option', opt);
end
Imap = reshape(uint8(COL(C,:)),h,w,3);

function ImHarris = encrust_harris(ImRGB, HC)
% encrust the Harris corners
% ImRGB : image source
% HC    : Harris corners

for r = 1:size(HC,1)
    ImRGB( HC(r,1)-2:HC(r,1)+2 , HC(r,2)-2           , 1) = 180;
    ImRGB( HC(r,1)-2:HC(r,1)+2 , HC(r,2)+2           , 1) = 180;
    ImRGB( HC(r,1)-2           , HC(r,2)-2:HC(r,2)+2 , 1) = 180;
    ImRGB( HC(r,1)+2           , HC(r,2)-2:HC(r,2)+2 , 1) = 180;

    ImRGB( HC(r,1)-2:HC(r,1)+2 , HC(r,2)-2           , 2) = 0;
    ImRGB( HC(r,1)-2:HC(r,1)+2 , HC(r,2)+2           , 2) = 0;
    ImRGB( HC(r,1)-2           , HC(r,2)-2:HC(r,2)+2 , 2) = 0;
    ImRGB( HC(r,1)+2           , HC(r,2)-2:HC(r,2)+2 , 2) = 0;

    ImRGB( HC(r,1)-2:HC(r,1)+2 , HC(r,2)-2           , 3) = 0;
    ImRGB( HC(r,1)-2:HC(r,1)+2 , HC(r,2)+2           , 3) = 0;
    ImRGB( HC(r,1)-2           , HC(r,2)-2:HC(r,2)+2 , 3) = 0;
    ImRGB( HC(r,1)+2           , HC(r,2)-2:HC(r,2)+2 , 3) = 0;

    ImRGB( HC(r,1)             , HC(r,2)             , 1) = 255;
    ImRGB( HC(r,1)             , HC(r,2)             , 2) = 0;
    ImRGB( HC(r,1)             , HC(r,2)             , 3) = 0;
end

ImHarris = ImRGB;

function disp_har(dbg, Ih, Hc)
% display the Harris corners
% Ih : Harris image
% Hc : Harris corners list
tith = sprintf('Harris corners - n = %d', size(Hc,1));

f5 = figure; set(gcf,'Color',[0.2,0.2,0.2]);
if dbg, set(f5,'Position',[2570 556 635 445]);
else    set(f5,'Position',[1680 0 1280 920]); end;

image(Ih); title(tith,'Color','w'); axis off; axis image;

function disp_abs_grad(dbg, Iax, Iay)
% display the absolute gradient images
IAX = uint8(255 * Iax / max(Iax(:)) );
IAY = uint8(255 * Iay / max(Iay(:)) );

f6 = figure(6); set(gcf,'Color',[0.2,0.2,0.2]);
if dbg, set(f6,'Position',[3210 10 630 440]);
else    set(f6,'Position',[1680 0 1280 920]); end;

subplot(2,1,1); image(repmat(IAX,[1 1 3]));
title('Iax','Color','w'); axis off; axis image;
subplot(2,1,2); image(repmat(IAY,[1 1 3]));
title('Iay','Color','w'); axis off; axis image;

function [NXY, AXY] = disp_norm_angle(dbg, Nxy, Axy)
% display the norm and the angle of the oriented x,y gradients
NXY = uint8( 255 * Nxy / max(Nxy(:)) );
AXY = uint8( 255 * (Axy + pi) / (2*pi) );

f7 = figure(7); set(gcf,'Color',[0.2,0.2,0.2]);
if dbg, set(f7,'Position',[3210 10 630 440]);
else    set(f7,'Position',[1680 0 1280 920]); end;

subplot(2,1,1); image(repmat(NXY,[1 1 3]));
axis image; axis off; title('Norm','Color','w');
subplot(2,1,2); image(repmat(AXY,[1 1 3]));
axis image; axis off; title('Angle','Color','w');

function [histw, histv] = histwv(v, w, vmin, vmax, bins)
% Inputs:
% v - values
% w - weights
% vmin - minimum value
% vmax - maximum value
% bins - number of bins (inclusive)

% Outputs:
% histw - weighted histogram
% histv (optional) - histogram of values

delta = (vmax-vmin)/(bins-1);
subs  = round((v-vmin)/delta)+1;

histv = accumarray(subs,1,[bins,1]);
histw = accumarray(subs,w,[bins,1]);

function disp_hist(dbg, histw, histv)
% display the histograms of the orientations
f8 = figure(8); set(gcf,'Color',[0.2,0.2,0.2]);
set(f8,'Position',[395 530 640 415]);
if dbg, set(f8,'Position',[3210 570 630 440]);
else    set(f8,'Position',[1680 0 1280 920]); end;

subplot(2,1,1); bar(histw); title('weighted histogram','Color','w');
subplot(2,1,2); bar(histv); title('normalized histogram','Color','w');

function disp_can(dbg, Ie)
% display the Canny edges
f9 = figure(9); set(gcf,'Color',[0.2,0.2,0.2]);
if dbg, set(f9,'Position',[3210 570 630 440]);
else    set(f9,'Position',[1680 0 1280 920]); end;
image(repmat(Ie,[1 1 3]));title('Canny','Color','w'); axis off; axis image;
