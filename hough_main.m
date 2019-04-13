% Hough main function - 04/16/2019
% Michel Dhome - Jean-Marc Berthommé
%
% - 04/11/2019:
%   . first script to handle Hough transform
% - 04/16/2019:
%   . deep work to manage binary images at multi-(theta,d) resolution

function hough_main
dbg = 0;      % debug flag
tres = 40;    % theta resolution - 19 - 180 between 180°
dres = 50;    % delta resolution - 21 - 100

% I. Build a random image
% h = 6; w = 8; % 6x8
h = input('image height ? ');
w = input('image width  ? ');
t = sprintf('Binary image: %dx%d', h,w);

Irnd = rand(h,w); % random matrix/image
I = Irnd <= 0.5;  % binary matrix/image

% I_ = zeros(h,w); I_(6,8) = 1;
% I = logical(I_);

% II. Display it
f1 = figure(1); set(gcf,'Color',[0.2,0.2,0.2]);
set(f1,'Position', [1316 532 605 442]);
image(repmat(I,[1 1 3])); axis image; axis off;
title(t,'color','w');

% White pixels coordinates
[X,Y] = meshgrid(1:w,1:h);
xx = X(I); yy = Y(I);
% [yy(1:10) xx(1:10)]

% Dual space [\theta, d] building
dt = -10;             % theta delta - [°] - 10°
tm =   0;             % theta   min - [°]
tM = 180;             % theta   max - [°] - 360 -> 180
tmin = tm + dt;       % real theta min - [°]
tmax = tM + dt;       % real theta max - [°]

dd = 0;                    % dist delta  - [px]
dmax = (h^2+w^2)^0.5 + dd; % dist max    - [px]
dmin =         -dmax + dd; % dist min    - [px]

% tres = 5;          % /!\ theta resolution - 19 - 180
% dres = 5;          % /!\ delta resolution - 21 - 100

t = linspace(tmin, tmax, tres)'; % angles vector - [°]
theta = t/180*pi;                 % thetas vector - [rad]

tstep = ( tmax-tmin ) / (tres - 1);
dstep = ( dmax-dmin ) / (dres - 1);

accu = zeros(dres, tres);
f2 = figure(2); set(f2,'Position', [6 9 1295 965]);
axis([tmin-tstep/2, tmax+tstep/2, -dmax-dstep/2, dmax+dstep/2]);
axis square;
disp_grid(tmin,tmax,tstep,tres, -dmax,dmax,dstep,dres, ':k');

nb = size(xx,1);
dedges = linspace(dmin-dstep/2, dmax+dstep/2, dres+1);

Iacc = uint8(255*accu);
f3 = figure(3); set(f3,'Position', [1311 11 605 434]);
image(repmat(Iacc, [1 1 3])); axis image;

Xd = (1:tres)'; % theta x coordinates

%**************** MAIN METHOD ***************
for i=1:nb
   y=yy(i); x=xx(i);               %     ->  ->
   d = -x*sin(theta)+y*cos(theta); % d = OM . e_\theta
   
   figure(2); hold on; plot(t, d); hold off;
   if (dbg==1), hold on; plot(t, d, 'r+'); hold off; end;
   
   [~, dBin] = histc(d, dedges);
   Yd = (dres+1)-dBin;        % theta y coordinates
   Ind = Yd + (Xd-1)* dres;   % Index of all the [y x] sinusoïde pixels
   accu(Ind) = accu(Ind) + 1; % increment the accumulation 
   
   % renormalizd accumulation image 
   Iacc = uint8(255*(accu/max(accu(:))));
   figure(3); image(repmat(Iacc, [1 1 3])); axis image;
end

pause; clear all; close all;

function disp_grid(xm, xM, dx, nx, ym, yM, dy, ny, color)
% display a grid
h = linspace(xm-dx/2, xM+dx/2, 20);
v = linspace(ym-dy/2, yM+dy/2, 20);

X = xm + dx * (0:nx)' - dx/2;
Y = ym + dy * (0:ny)' - dy/2;

hold on;
for j = 1:nx+1     % vertical lines
    x = X(j);
    plot(x*ones(1,20), v, color);
end

for i = 1:ny+1     % horizontal lines
    y = Y(i);
    plot(h, y*ones(1,20), color);
end
hold off;

