    function plaidData = plaidDataPerFrame(deg_per_px, ifi, pa, x, y, u, v)
    sf =.1; % cycles per deg
    sf = (deg_per_px).*sf; %1/px_per_deg %convert from deg to 
    af = 2*pi*sf;
    Vx = (1/deg_per_px).*u; %(1/mmPerPixel)*(tand(u/2)*2*viewingDistance); % (dg/s) -> px/s
    Vy = (1/deg_per_px).*v;%1/mmPerPixel)*(tand(v/2)*2*viewingDistance); 

    pa.dtx = pa.dtx + ((2 * pi * (Vx * (1/1000) * sf)) / ifi); % px/s * s/ms * cyc/px * ms/frame = cyc/frame
    pa.dty = pa.dty + ((2 * pi * (Vy * (1/1000) * sf)) / ifi);
    xGrating = (.5 + .5 * sin(af * x - pa.dtx)); % vertically oriented grating
    yGrating = (.5 + .5 * cos(af * y - pa.dty)); % horiz grating
    plaidData = 255.*((xGrating + yGrating)/2);
    plaidData = repmat(plaidData,[ 1 1 3 ]);
    %x = single(plaidData);
    plaidData = permute(plaidData,[ 3 2 1 ]);
    
    %{
    figure; imagesc(255.*((xGrating + yGrating)/2))
    axis equal
    colormap gray
    shg
    keyboard
%}
end