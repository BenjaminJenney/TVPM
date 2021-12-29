% FangFang

% This function takes the number of bins, the boundaries and the data as
%inputs and generates binned data

% GroupedData format:
% 1st row: stimulus intensity
% 2nd row: number of total trials tested at that level of intensity
% 3rd row: number of trials judging the C is to the right of the S
function GroupedData = BinData(numBins,bounds,data)
    %calculate the mid value of each interval
    scalePoints = linspace(bounds(1),bounds(2),numBins+1);
    Bins = (scalePoints(2:end)+scalePoints(1:end-1))./2;
    %tol is the distance between the interval boundaries and the mid value
    %of that bin
    tol = diff(Bins(1:2))/2;
    
    %initialize the matrix GroupedData
    GroupedData = zeros(3,length(Bins));
    %the first row is the stimulus intensity (the mid value of each interval)
    GroupedData(1,:) = Bins;
    for l = 1:length(Bins)
        %find the indices that correspond to the data falling into the bin
        indices_fallingOnTheBin = find(abs(data(1,:)-Bins(l))<=tol);
        %sum the number of those data
        GroupedData(2,l) = length(indices_fallingOnTheBin);
        %Given those trials, how many are "right" responses
        GroupedData(3,l) = sum(data(2,indices_fallingOnTheBin));
    end     
    %if no trial falls on a bin, delete that column so when plotting it,
    %it wouldn't make any error
    GroupedData(:,GroupedData(2,:)==0) = [];
end