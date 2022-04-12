# Optic-Flow-Stereo

3D Pong Experiment

To run in debug mode – without HMD connected, in SetupDisplayforSDK2 on line 7, set:

ds.oculusConnected = 0; % Is the HMD connected

For production mode, in RunDK2Experiment make sure to comment out:

Screen('Preference', 'SkipSyncTests', 1);

