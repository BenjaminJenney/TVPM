# TVPM

To run expirements NB: added flag to stream to headset or screen:

Params: (subjectInitials,condition,monocularFlag,trainingFlag,mask, streamToHMDFlag = 0 to stream to monitor, 1 to stream to HMD)

Reccomended inputs for running various experiments change last param to 0 or 1 depending if you want to stream to HMD or monitor:

Training real eye movement: RunTVPMExperiment('test', 'real', 0, 0, 0, 0)  

Training simulated eye movement: RunTVPMExperiment('test', 'simulated', 0, 0, 0, 0)

TVPM static disparity: RunTVPMExperiment('test', 'tvpmsd', 0, 0, 0, 0)

TVPM changing disparity: RunTVPMExperiment('test', 'tvpmcd', 0, 0, 0, 0)

Of note: you can control stimuli duration by changing pa.stimulusDuration_sec in SetupParameters.
