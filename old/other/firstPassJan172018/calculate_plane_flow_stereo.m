function [velocity_field] = calculate_plane_flow_stereo(XYZ, centers_m, translate, rotate, focalLength)

% CSB RENAME everything here so this is clear and consistent w
% make_dot_plane_stereo... !!! then test old script and see format of
% output and mimic that here as well.

% The world is a plane, the observer is translating and rotating past it.
% The observer is always at the origin, and is looking down the Z axis.
% This program takes a set of gabors and calculates the velocity of each
% one based on given translational and rotational components (for phase heading).
%
% XYZ is the XYZ position (in m) of the locations where we are measuring
% the speed. It represents the 3-D coordinates of the dots.
%
% centers_m is the x,y positions of the patches (dots) on the screen (in m).
%
% translate tells you the observer's translation speed (Vx, Vy, Vz) in
% m/sec.
%
% rotate tells you the observer's rotation (Wx, Wy, Wz) in rad/s.
%
% focalLength is how far the observer is the from the screen in meters.
%
% velocity_field is a matrix (2x[number of locations]) reporting
% the x and y velocity in deg/s.
%
% CSB
% updated Jan 19, 2018

if nargin==0
    [unused1,centers_m,unused2,XYZ] = make_dot_plane_stereo2(.1, 12.5, [40 32], [-20 -3 20 3],.57, 0.06285);
    translate = [0.14 0 1.9];                                               % the observer is walking forward at 1.9 m/s, a brisk walk
    rotate = [0 0 0];                                                       % no rotation
    focalLength = .57;
end 

% initialize empty matrix to hold velocities
velocity_field1 = zeros(size(centers_m));                                 % number of patches, X and Y vel

%solve for the inverse depth at the X,Y,Z position that projects to the
%desired x,y position, then use that to compute the flow vector at each
%location (same Z value for all of the gabors).
for i=1:length(XYZ)
        inv_depth = 1/XYZ(3,i);
        
        % All these equations come from Heeger and Jepson's chapter
        x = centers_m(1,i);
        y = centers_m(2,i);
        f = focalLength;
        Pxy = inv_depth;
        
        A = [-f 0 x; 0 -f y];
        B = [(x*y)/f -(f+x^2/f) y; f+y^2/f -(x*y)/f -x];
        
        velocity_field1(1:2,i) = Pxy*A*translate'+B*rotate';
end

%gonna convert velocities in m/s to deg/sec, because that's what we actually
%care about.
velocity_field = atand(velocity_field1/focalLength);

% to plot flow field, to debug:
%{
figure
quiver(centers_m(1,:),centers_m(2,:),velocity_field(1,:),velocity_field(2,:)) % m
figure
quiver(atand(centers_m(1,:)/focalLength),atand(centers_m(2,:)/focalLength),velocity_field(1,:),velocity_field(2,:)) % deg
%}

end