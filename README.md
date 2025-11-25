# linear-encoder-evaluation
This repository contains the MATLAB script used for the evaluation of the linear encoder developed during my TDK work.
The MATLAB script used for the measurement evaluation was developed with assistance from ChatGPT 5.1 Thinking, as documented in the agents.md configuration file within the project repository. After that, manual adjustments have been made to the script, which contains comments explaining what to do to enable certain functionalities.

eval_11_raw_100_hz.csv contains data recorded during evaluation. This file is handled by the MATLAB script. The firsc column is the measurement id. As samples were taken at a rate of 100 Hz, there are 100 rows per second. The second and third columns contain raw data sent by the sensor, while the last column contains ADC measurements of the photodiode signal.

raw20_cropped.csv contains data recorded during initial testing. This file does not contain interferometer data. The ladder was pulled quickly by hand for this recording, causing high-frequency sinusoidal waves.
