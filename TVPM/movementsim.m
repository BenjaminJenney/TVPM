%% computes heading and rotation velocity over time
% for fixating observer, with parameters from REM stereo pilot experiment
% charlie burlingham, june 29, 2021

theta0 = [5,10,15]; %deg. initial heading

for ii = 1:length(theta0)
    
    ly = 15; %m
    lx = ly.*tand(theta0(ii)); %m
    lz = sqrt(lx.^2 +ly.^2); %m
    s = 5; %m/s
    t = 0:.1:2; %s
    
    theta = asind(lx./(lz-s.*t)); % deg. heading over time
    rotVel = gradient(theta); % deg/s. rotation velocity over time
    
    subplot(2,1,1);
    hold on;
    plot(t,theta);
    ylabel('Heading (deg)')
    
    subplot(2,1,2);
    hold on;
    plot(t,rotVel);
    ylabel('Rotation velocity (deg/s)')
    xlabel('Time (s)')
end

legend({'heading = 5 deg', '10 deg' '15 deg'})