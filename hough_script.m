% Hough function - 04/18/2019
% Jean-Marc Berthommé
%
% - 04/11/2019:
%   . first draft to handle the Hough transform
% - 04/13/2019:
%   . deep work to manage binary images at multi-(theta,d) resolution
%   . code cleaning & factorization

function hough_script
dbg = 1;   % debug flag
tres = 19; % theta resolution: [19,180] for ~180°
dres = 21; % delta resolution: [21,150]

% Build & display a binary image
h = 6; w = 8; % 6x8 - 16x32
% h = input('image height ? '); w = input('image width  ? ');
Irnd = rand(h,w); % random matrix/image
I = Irnd < 0.5;   % binary matrix/image
% I_ = zeros(h,w); I_(6,8) = 1; I = logical(I_);
% I_ = ones(h,w);  I_(6,8) = 0; I = logical(I_);
fig = disp_binary_img(I);

% *** CALL THE MAIN METHOD ***
[accu, t, d] = hough(I, tres, dres, dbg, fig);

% Find the best line
[accu_max, id_max] = max(accu(:));
[T,D] = meshgrid(1:tres,1:dres);
tb = T(id_max); db = D(id_max);

fprintf('Best line accumulated %d pixels', accu_max);
fprintf(' for theta = %0.3f° and d = %0.3f px.\n', t(tb), d(db));
fprintf('This matches to pixel [%d, %d].\n', tb, db);

figure(fig+2); hold on; plot(tb, db, 'r.'); hold off;
figure(fig+1); hold on; plot(t(tb), d(db), 'r.'); hold off;

pause; clear all; close all;

function f1 = disp_binary_img(I)
% display a binary image
[h, w] = size(I);
t1 = sprintf('Binary image: %dx%d', h,w);
f1 = figure(1); set(gcf,'Color',[0.2,0.2,0.2]);
set(f1,'Position', [1316 532 605 442]);
image(repmat(I,[1 1 3])); axis image; axis off; title(t1,'color','w');