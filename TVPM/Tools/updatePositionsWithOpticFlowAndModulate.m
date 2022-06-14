function [plaidPosX, plaidPosY] = updatePositionsWithOpticFlowAndModulate(numPlanes, currentTrialNumber, frameCount, sineCycle_m, plaidPosX, plaidPosY, vXw, vYw, initialPlaidPosX_m, initialPlaidPosY_m)
    for i = numPlanes
        plaidPosX{i} = plaidPosX{i} + squeeze(vXw(currentTrialNumber, frameCount, i, :)); % update previous plaid positions with the displacements generated from preprocessing (Optic Flow)
        plaidPosY{i} = plaidPosY{i} + squeeze(vYw(currentTrialNumber, frameCount, i, :));

       % ATTN: MAKE MODULATION MODULAR
       % attempt at modulation
        displacementsX{i} = initialPlaidPosX_m{i}-plaidPosX{i};
        displacementsY{i} = initialPlaidPosY_m{i}-plaidPosY{i};

        for jj = 1:length(plaidPosX{i})
            if displacementsX{i}(jj) >= sineCycle_m(i)
                plaidPosX{i}(jj) = initialPlaidPosX_m{i}(jj) + -sign(displacementsX{i}(jj))*(abs(displacementsX{i}(jj))-sineCycle_m(i)); % this modulates to the phase beyond the threshold if it crosses the threshold, preventing jumps
            end

            if displacementsY{i}(jj) >= sineCycle_m(i)
                plaidPosY{i}(jj) = initialPlaidPosY_m{i}(jj) + -sign(displacementsY{i}(jj))*(abs(displacementsY{i}(jj))-sineCycle_m(i));
            end
        end


    end
end