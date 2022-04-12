%% Practice Bias Function
 
gazeRate = [0 1 2 3 4 5 6 7 8 9];
activeBias = [0.05 0.135 0.21 0.23 0.25 0.254 0.256 0.262 0.293 0.336];
simulatedBias = [0.35 0.535 0.81 1.33 1 1.54 1.65 1.82 2 2.3];
simulatedBias([2,3,4,5,6,7,8,10]) = simulatedBias([2,3,4,5,6,7,8,10])+0.336
simulatedBias([1,8,9]) = simulatedBias([1, 8,9])+0.256;
plot(gazeRate,activeBias, '-o')
hold on
plot(gazeRate, simulatedBias, '-o')
hold off
ylim([-3,3])
ylabel('Bias^{\circ}', 'FontSize',20,'FontWeight','bold')
xlabel('Gaze rotation rate ^{\circ}/sec', 'FontSize',20,'FontWeight','bold')
title('Sample Bias Function', 'FontSize', 28, 'FontWeight', 'bold')
legend('Location','southeast')
legend('Active Head Pursuit', 'Simulated Head Pursuit', 'FontSize',16)
matlab2tikz('WorkingBiasFunction.tex')
 
%% Practice Bias Function -- Both Conditions Failed
 
gazeRate = [0 1 2 3 4 5 6 7 8 9];
activeBias = [0.35 0.535 0.81 1.33 1 1.54 1.65 1.82 2 2.3];
simulatedBias = [0.42 0.57 0.79 1.22 1.3 1.57 1.78 1.98 2.12 2.56];
plot(gazeRate,activeBias, '-o')
hold on
plot(gazeRate, simulatedBias, '-o')
hold off
ylim([-3,3])
ylabel('Bias^{\circ}', 'FontSize',20,'FontWeight','bold')
xlabel('Gaze rotation rate ^{\circ}/sec', 'FontSize',20,'FontWeight','bold')
title('Sample Bias Function -- Both Conditions Biased', 'FontSize', 28, 'FontWeight', 'bold')
legend('Location','southeast')
legend('Active Head Pursuit', 'Simulated Head Pursuit', 'FontSize',16)
matlab2tikz('NotWorkingBiasFunction.tex')
 
%% Practice Bias Function -- No bias
 
gazeRate = [0 1 2 3 4 5 6 7 8 9];
activeBias = [0.05 0.135 0.21 0.23 0.25 0.254 0.256 0.262 0.293 0.336];
simulatedBias = activeBias;
simulatedBias([2,3,4,5,6,7,8]) = simulatedBias([2,3,4,5,6,7,8])+0.016;
simulatedBias([1,8,9]) = simulatedBias([1, 8,9])+0.056;
simulatedBias([10]) = simulatedBias([10]) - 0.042;
plot(gazeRate,activeBias, '-o')
hold on
plot(gazeRate, simulatedBias, '-o')
hold off
ylim([-3,3])
ylabel('Bias^{\circ}', 'FontSize',20,'FontWeight','bold')
xlabel('Gaze rotation rate ^{\circ}/sec', 'FontSize',20,'FontWeight','bold')
title('Sample Bias Function -- No Bias', 'FontSize', 28, 'FontWeight', 'bold')
legend('Location','southeast')
legend('Active Head Pursuit', 'Simulated Head Pursuit', 'FontSize',16)
matlab2tikz('NoBiasFunction.tex')

 