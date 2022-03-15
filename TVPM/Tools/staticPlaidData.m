function plaidData = staticPlaidData(ds, sf, x, y)
    white = 1.0;
    grey  = white/2;
    inc   = white - grey;
    
    %sf = 2; % cycles per deg .1 * (the ratio between the largest dimension of disk texture and the largest dimension of the plane texture)
    sf = (ds.deg_per_px).*sf; %1/px_per_deg %convert from deg to 
    af = 2*pi*sf;
    
    xGrating = (grey + inc * sin(af * x)); % vertically oriented grating
    yGrating = (grey + inc * cos(af * y)); % horiz grating
    plaidData = 255.*((xGrating + yGrating)/2);
    plaidData = repmat(plaidData,[ 1 1 3 ]);
    plaidData = permute(uint8(plaidData),[ 3 2 1 ]);
    plaidData(4,:,:) = 255;
    %{
    figure; imagesc(255.*((xGrating + yGrating)/2))
    axis equal
    colormap gray
    shg
    keyboard
%}
end