# Orbital-Propagator
Orbital Propagator made for AE502 using universal anomaly method.

# Quick, informal summary of methodology for development and testing
## Universal Propagator Testing
I tested my universal propagator against the universal propagator in the python package “pykep”. In the test, I propagated a solar system consisting of the two bodies, Earth and Mars, for approximately a year and found that there was an error < 10e-4 for the entire run.

## Izzio Lambert Solver Testing
Using Izzio’s formulation, I developed my own script for solving the lambert problem. Similarly to the universal propagator, I tested my code against the lambert solver in “pykep”. Both my code and the lambert solver from pykep were executed against over 300,000 custom Earth to Mars transfer scenarios developed using ephemeris data from JPL’s Horizons tool and my recently written and tested universal propagator script. Out of ~320,000 runs, only ~11 instances were found where the difference between the pykep output and my output ever exceeded 1e-5.



## Porkchop Plots
![Earth to Borisov Rendevous](https://user-images.githubusercontent.com/92574647/221439369-07635d4e-b159-410d-afa0-6e44f87a2b82.png)
![Earth to Oumouamoua Rendevous](https://user-images.githubusercontent.com/92574647/221439371-51456f0b-f4dc-4cf0-bb50-a7b40c4be592.png)
![Earth to Oumouamoua Fly-by](https://user-images.githubusercontent.com/92574647/221439370-bab4edfd-6ccd-4a95-b552-7a22fb1ce0b8.png)
![Earth to Borisov Fly-by](https://user-images.githubusercontent.com/92574647/221439372-70010d54-6f9c-49f0-af98-20e963eb3b00.png)

![Orbital Element Table](https://user-images.githubusercontent.com/92574647/221440260-384993c7-6e69-401d-899b-4593d26ac5d4.PNG)
