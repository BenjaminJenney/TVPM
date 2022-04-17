# TVPM

To run expirements:
Training real eye movement: RunTVPMExperiment('test', 'real', 0, 0, 0)
Training simulated eye movement: RunTVPMExperiment('test', 'simulated', 0, 0, 0)
TVPM static disparity: RunTVPMExperiment('test', 'tvpmsd', 0, 0, 0)
TVPM changing disparity: RunTVPMExperiment('test', 'tvpmcd', 0, 0, 0)

Of note: you can control stimuli duration by changing pa.stimulusDuration in SetupParameters.
