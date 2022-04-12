function [dots2D_L, dots2D_R, dots_3D] = make_dot_plane_stereo(density, dist, dims, exclude, focalLength)
% Make a plane of random dots.
%
% density tells you the dot density in dots/deg^2 for surfaces.
%        Density from Warren & Hannon (1990) works out to about .22
%        dots/deg^2.
%
% dist tells you how far the observer starts from a dot surface in m.
%
% dims tells you the [X,Y] dimensions of the dot field on the screen in
%       deg. in deg.
%
% exclude is a rectangle [X1, Y1, X2, Y2] where no dots will be placed,
%       lets you do things like block out the FOE. Can be set to 0 to
%       display the whole world. in deg
%
% dots2D_L is a three row vector containing the [X;Y;Z] positions of the
%       dots, in m, for the left eye
%
% dots2D_R is a three row vector containing the [X;Y;Z] positions of the
%       dots, in m, for the right eye
%
% dots2D_L is a two row vector containing the [X;Y] positions of the dots
%       on the screen, in m, for the left eye
%
% dots2D_R is a two row vector containing the [X;Y] positions of the dots
%       on the screen, in m, for the right eye

if nargin==0 % for testing & debug
    density = .1; % dots/deg^2
    dist = 12.5; % m
    dims = [40 32]; % deg
    exclude = [-20 -3 20 3]; % deg
    focalLength = .57; % m
    IPD = 0.06285; % meters. This is the US average btw men and woman
end

% First take the dimensions of the screen (same as the dims of the plane, but doesn't have to be) and the
% "exclude" region in degrees and compute a plane (with rect cutout) composed of randomly positioned dots. 
nDots = ceil(prod(dims)*density);
dots_deg_temp = [(rand(1,nDots)-.5)*dims(1); (rand(1,nDots)-.5)*dims(2)];
legal = ones(1,nDots);
if exclude
    for i=1:nDots
        if IsInRect(dots_deg_temp(1,i),dots_deg_temp(2,i),exclude)
            legal(i) = 0;
        end
    end
end

eyeShiftAmount = IPD/2; % m. 
f = focalLength; % distance of actual observer's eye to the display (HMD screen or CRT monitor)
Z = dist; % distance of simulated observer / camera to the plane

dots_2D = f.*tand(dots_deg_temp(:,legal==1));
dots_3D(1:2,:) = Z .* tand(dots_2D); % now in m. Convert only the ones that aren't in exclusion region
dots_3D(3,:) = Z;

X = dots_3D(1,:);
Y = dots_3D(2,:);

x_R = ((f.*X)/Z) + eyeShiftAmount;
x_L = ((f.*X)/Z) - eyeShiftAmount;

%{
X = dots_3D(1,:);
Y = dots_3D(2,:);

x_R = ((f*X)/Z) + eyeShiftAmount;
x_L = ((f*X)/Z) - eyeShiftAmount;

y = ((f*Y)/Z);

dots2D_R = [x_R ; y]; % deg
dots2D_L = [x_L ; y]; % deg
%}




dots2D_R = f*tand(dots2D_R); % convert to m
dots2D_L = f*tand(dots2D_L); % convert to m

%}
%{
IPD = 0.06285; % meters. This is the US average btw men and woman
eyeShiftAmount = IPD/2; % m. 

f = focalLength; % distance of actual observer's eye to the display (HMD screen or CRT monitor)
Z = dist; % distance of simulated observer / camera to the plane

% convert deg to m in simulated image plane space
dots_2D = f.*tand(dots_deg_temp(:,legal==1)); % now in m. Convert only the ones that aren't in exclusion region

dots_3D = (dots_2D.*Z)./f; % m. 
dots_3D(3,:) = Z; % add last dimension, which is constant, because the plane is frontoparallel to camera

X = dots_3D(1,:);
Y = dots_3D(2,:);

x_R = ((f*X)/Z) + eyeShiftAmount;
x_L = ((f*X)/Z) - eyeShiftAmount;

y = ((f*Y)/Z);

dots2D_R = [x_R ; y]; % m
dots2D_L = [x_L ; y]; % m
%}

%{
IPD = 4;
eyeShiftAmount = IPD/2; % degrees % This should be half the IPD- as we measure physically for each subject when they fixate at a very far away point

dots3D_L(1,:) = dots_deg_temp(1,:) - eyeShiftAmount; dots3D_L(2,:) = dots_deg_temp(2,:); 
dots3D_R(1,:) = dots_deg_temp(1,:) + eyeShiftAmount; dots3D_R(2,:) = dots_deg_temp(2,:);


dots3D_L = dots3D_L(:,legal==1);
dots2D_L(1:2,:) = dist*tand(dots3D_L); % deg to m, in the virtual 3-D space, hence why using the distance from simulated cam to plane
dots2D_L(3,:) = dist;

dots3D_R = dots3D_R(:,legal==1);
dots2D_R(1:2,:) = dist*tand(dots3D_R); % deg to m
dots2D_R(3,:) = dist;

% Convert all image plane dot positions to m, in real world.
dots2D_L = focalLength*tand(dots2D_L); % deg to m, but in real 3-D space, hence why using the focal length ? view dist from observer eye to display. 
dots2D_R = focalLength*tand(dots2D_R);
%}