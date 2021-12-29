function [movieData, numFrames] = moviePlaidData(ds, pa, x, y, u, v)
    white = 1.0;
    grey  = white/2;
    inc   = white - grey;
    
    numFrames = ceil((1/ds.ifi) * pa.targetMotionDuration);
    sf =.1; % cycles per deg
    sf = (ds.deg_per_px).*sf; %1/px_per_deg %convert from deg to 
    af = 2*pi*sf;
    Vx = (1/ds.deg_per_px).*u; %(1/mmPerPixel)*(tand(u/2)*2*viewingDistance); % (dg/s) -> px/s
    Vy = (1/ds.deg_per_px).*v;%1/mmPerPixel)*(tand(v/2)*2*viewingDistance); 
    
    movieData = {}; dtx = 0; dty = 0;
    for i = 1:numFrames 
        dtx = dtx + ((2 * pi * (Vx * (1/1000) * sf)) / ds.ifi); % px/s * s/ms * cyc/px * ms/frame = cyc/frame
        dty = dty + ((2 * pi * (Vy * (1/1000) * sf)) / ds.ifi);
        xGrating = (grey + inc * sin(af * x - dtx)); % vertically oriented grating
        yGrating = (grey + inc * cos(af * y - dty)); % horiz grating
        plaidData = 255.*((xGrating + yGrating)/2);
        plaidData = repmat(plaidData,[ 1 1 3 ]);
        plaidData = permute(uint8(plaidData),[ 3 2 1 ]);
        movieData{i} = plaidData;
    end
    
    
    %{
    figure; imagesc(255.*((xGrating + yGrating)/2))
    axis equal
    colormap gray
    shg
    keyboard
%}
end