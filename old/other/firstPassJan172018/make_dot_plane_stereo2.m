function [dots2D_L, dots2D_R, dots3D_L, dots3D_R] = make_dot_plane_stereo2(density, dist, dims, exclude, focalLength ,IPD)
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
    IPD = 0.06285; % meters. measure this physically for each observer
end

eyeShiftAmount = IPD/2; % m.
f = focalLength; % distance of actual observer's eye to the display (HMD screen or CRT monitor)
Z = dist; % distance of simulated observer / camera to the plane

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

dots2D_L(1,:) = ( tand(dots_deg_temp(1,:)).*f ) - eyeShiftAmount; % x coordinates (translated) 
dots2D_L(2,:) = ( tand(dots_deg_temp(2,:)).*f ); % y coordinates (not translated)
dots2D_R(1,:) = ( tand(dots_deg_temp(1,:)).*f ) + eyeShiftAmount; % x coordinates
dots2D_R(2,:) = ( tand(dots_deg_temp(2,:)).*f ); % y coordinates 

dots3D_L(1,:) = ( tand(dots_deg_temp(1,:)).* (f + Z) ) - eyeShiftAmount; % x coordinates (translated) 
dots3D_L(2,:) = ( tand(dots_deg_temp(2,:)).* (f + Z) );  % y coordinates (not translated)
dots3D_R(1,:) = ( tand(dots_deg_temp(1,:)).* (f + Z) ) + eyeShiftAmount;  % x coordinates
dots3D_R(2,:) = ( tand(dots_deg_temp(2,:)).* (f + Z) ); % y coordinates 

dots3D_L(3,:) = Z;
dots3D_R(3,:) = Z;
