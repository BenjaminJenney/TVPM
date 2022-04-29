# TVPM

To run experiments:

TODO: What argument should be added for the planes variable, 3? So:

RunTVPMExperiment('test', 'real', 0, 0, 0, 3)?

Training real eye movement: RunTVPMExperiment('test', 'real', 0, 0, 0)

Training simulated eye movement: RunTVPMExperiment('test', 'simulated', 0, 0, 0)

TVPM static disparity: RunTVPMExperiment('test', 'tvpmsd', 0, 0, 0)

TVPM changing disparity: RunTVPMExperiment('test', 'tvpmcd', 0, 0, 0)

Of note: you can control stimulus duration by changing pa.stimulusDuration_sec in SetupParameters.
