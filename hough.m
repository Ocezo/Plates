% Hough function - 04/20/2019
% Jean-Marc Berthommé
%
% - 04/11/2019:
%   . first draft to handle the Hough transform
% - 04/13/2019:
%   . deep work to manage binary images at multi-(theta,d) resolution
%   . code cleaning & factorization
% - 04/18/2019:
%   . new protection against pixels outside the accumulation grid (~> Ok)
%   . modularity added to call a function in the plate_recognition project
%   . random subsampling of the detections to fasten the "accu" estimation
% - 04/20/2019:
%   . "accu" estimation is now ~x300 faster ("for" ~> matrix calculations)

function [accu, t, d] = hough(I, tres, dres, max_nb, dbg, fig)
% Hough transform of a binary image at resolution [dres, tres]
% I      : binary input image
% tres   : theta    resolution - px (parameter dual space)
% dres   : distance resolution - px (parameter dual space)
% max_nb : targeted number of detections (~> random shrinkage)
% dbg    : debug flag
% fig    : current figure number

% Get the detections matching to the white pixels
[h,w] = size(I);
[X,Y] = meshgrid(1:w,1:h);
xx = X(I); yy = Y(I);

% Violent random subsampling to accelerate the "accu" estimation
nb = size(xx,1); % nb of white pixels i.e. nb of detections

if  nb > max_nb
    rnd_ = rand(nb,1);
    rnd = rnd_ < max_nb/nb;
    fprintf('%d white pixels are picked from %d ', sum(rnd), nb);
    fprintf('to accelerate the estimation of\nthe accumulation matrix.\n');
    
    % white pixels number shrinkage
    % xx2 = xx(rnd); yy2 = yy(rnd);
    xx = xx(rnd); yy = yy(rnd); % <- information loss!
end

% Check the "sieved" image
% I2 = zeros(h,w);
% for i=1:size(xx,1), I2(yy(i),xx(i))=1; end;
% figure; image(repmat(I2, [1 1 3])); axis image;
nb = size(xx,1); % new real nb of detections considered

% Build the parameter space [theta, d]
dt = -10;       % theta    delta - [°] - 10°
tm =   0;       % theta      min - [°]
tM = 180;       % theta      max - [°] - 360 -> 180 (pi-periodic)
tmin = tm + dt; % real theta min - [°]
tmax = tM + dt; % real theta max - [°]
% tres = 19;    % /!\ theta resolution -  [°] - 19/180
% tstep = ( tmax-tmin ) / (tres - 1);

t = linspace(tmin, tmax, tres)'; % angle theta vector -  [°]
theta = t/180*pi;                % angle theta vector - [rad]

dd = 0;             % dist delta    - [px]
dM = (h^2+w^2)^0.5; % dist max      - [px]
dmin = -dM + dd;    % real dist min - [px]
dmax =  dM + dd;    % real dist max - [px]
% dres = 21;        % /!\ delta resolution - [px] - 21/100
dstep = ( dmax-dmin ) / (dres - 1);

d = linspace(dmax, dmin, dres)'; % dist vector - [px]
dedges = linspace(dmin-dstep/2, dmax+dstep/2, dres+1);

% Accumulation matrix
accu = zeros(dres, tres);

% First displays
% disp_first_dual(tmin, tmax, tres, dmin, dmax, dres, fig);
t3 = sprintf('Dual space [theta, d]: %dx%d', dres, tres);
% disp_first_accu(accu, t3, fig);

% Xd = (1:tres)';  % theta x coordinates
% *** MAIN METHOD ***
% tic;
% for i=1:nb % TODO: hugely fasten this slow "for"! :$ ;D
%    y=yy(i); x=xx(i);                %     ->  ->
%    d_ = -x*sin(theta)+y*cos(theta); % d_ = OM . e_\theta
% 
%    [~, dBin] = histc(d_, dedges);   % make fall the points into the bins
%    Yd = (dres+1)-dBin;              % theta y coordinates
%    Ok = dBin ~= 0;                  % remove the pixels outside the grid
%    Ind = Yd(Ok) + (Xd(Ok)-1)* dres; % indexes of the sinusoïde pixels
%    accu(Ind) = accu(Ind) + 1;       % increment the accumulation matrix
%    
%    % display the dual space & the renormalized accumulation image
%    % disp_dual(t, d_, dbg, fig);
%    % disp_accu(accu, t3, fig);
% end
% toc;

tic;
% T = repmat(t', [nb, 1]);
D = -xx*sin(theta')+yy*cos(theta');
[~, Dbin] = histc(D, dedges);
OOk = Dbin ~= 0;
YYd = (dres+1)-Dbin;
XXd = repmat(1:tres, [nb 1]);
IND = YYd(OOk) + (XXd(OOk)-1)* dres;
% IND = reshape(IND, nb, tres);
ntot = dres*tres;
% N = histc(IND(:), (1:ntot) );
N = histc(IND, (1:ntot) );
accu2 = reshape(N, [dres, tres]);
toc;

% Check we have the same output
% accu_diff = accu - accu2;
% sum(accu_diff(:)) % = 0 ! ;D

% Final displays
disp_final_accu(accu2, t3, fig);
% disp_final_dual(tmin,tmax,tres, dmin,dmax,dres, fig, T, D, dbg);

% function disp_first_dual(tmin, tmax, tres, dmin, dmax, dres, fig)
% % first display of the dual space
% tstep = ( tmax-tmin ) / (tres - 1);
% dstep = ( dmax-dmin ) / (dres - 1);
% 
% f2 = figure(fig+1); set(f2,'Position', [6 9 1295 965]);
% axis([tmin-tstep/2,tmax+tstep/2, dmin-dstep/2,dmax+dstep/2]); axis square;
% disp_grid(tmin,tmax,tres,  dmin,dmax,dres, ':k');

function disp_grid(xm,xM,nx, ym,yM,ny, color)
% display a grid
n = 20; % line resolution

dx = (xM-xm) / (nx -1);
dy = (yM-ym) / (ny -1);
X = xm + dx * (0:nx)' - dx/2;
Y = ym + dy * (0:ny)' - dy/2;
h = linspace(xm-dx/2, xM+dx/2, n); % horizontal
v = linspace(ym-dy/2, yM+dy/2, n); % vertical

hold on;
for j = 1:nx+1 % vertical lines
    x = X(j);
    plot(x*ones(1,n), v, color);
end
for i = 1:ny+1 % horizontal lines
    y = Y(i);
    plot(h, y*ones(1,n), color);
end
hold off;

% function disp_dual(t, d_, dbg, fig)
% % display of the dual space
% figure(fig+1); hold on; plot(t, d_); hold off;
% if (dbg==1), hold on; plot(t, d_, 'm+'); hold off; end;

function disp_final_dual(tmin,tmax,tres, dmin,dmax,dres, fig, T, D, dbg)
% final display of the dual space
tstep = ( tmax-tmin ) / (tres - 1);
dstep = ( dmax-dmin ) / (dres - 1);

f2 = figure(fig+1); set(f2,'Position', [6 9 1295 965]);
axis([tmin-tstep/2,tmax+tstep/2, dmin-dstep/2,dmax+dstep/2]); axis square;
disp_grid(tmin,tmax,tres,  dmin,dmax,dres, ':k');

hold on; plot(T', D', 'b'); hold off;
if (dbg==1), hold on; plot(T', D', 'm+'); hold off; end;

% function disp_first_accu(accu, t3, fig)
% % first display of the accumulation matrix
% Iacc = uint8(255*accu);
% 
% f3 = figure(fig+2); set(gcf,'Color',[0.2,0.2,0.2]);
% set(f3,'Position', [1311 11 605 434]);
% image(repmat(Iacc, [1 1 3])); axis image; title(t3,'color','w');

% function disp_accu(accu, t3, fig)
% % display the accumulation matrix
% Iacc = uint8(255*(accu/max(accu(:))));
% 
% figure(fig+2); image(repmat(Iacc, [1 1 3])); axis image;
% title(t3, 'color','w');

function disp_final_accu(accu, t3, fig)
% display the final accumulation matrix
Iacc = uint8(255*(accu/max(accu(:))));

f3 = figure(fig+2); set(gcf,'Color',[0.2,0.2,0.2]);
set(f3,'Position', [1311 11 605 434]);
image(repmat(Iacc, [1 1 3])); axis image; title(t3,'color','w');

