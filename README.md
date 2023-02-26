# Orbital-Propagator
Orbital Propagator made for AE502 using universal anomaly method.
*IMPORTANT*
I don’t know how to work with branches on github right now and I’ve accidentally stored my code in the master branch and not the main branch.


# Quick, informal summary of methodology for development and testing
## Universal Propagator Testing
I tested my universal propagator against the universal propagator in the python package “pykep”. In the test, I propagated a solar system consisting of the two bodies, Earth and Mars, for approximately a year and found that there was an error < 10e-4 for the entire run.

## Izzio Lambert Solver Testing
Using Izzio’s formulation, I developed my own script for solving the lambert problem. Similarly to the universal propagator, I tested my code against the lambert solver in “pykep”. Both my code and the lambert solver from pykep were executed against over 300,000 custom Earth to Mars transfer scenarios developed using ephemeris data from JPL’s Horizons tool and my recently written and tested universal propagator script. Out of ~320,000 runs, only ~11 instances were found where the difference between the pykep output and my output ever exceeded 1e-5. Therefore, more debugging and testing is needed before this tool is deemed acceptable, but because the occurrances were so rare and most of the errors were less than 5% error with the lasgest being ~11%, I still used it in my porkchop plot creation.



## Porkchop Plots
![Earth to Borisov Rendevous](https://user-images.githubusercontent.com/92574647/221439369-07635d4e-b159-410d-afa0-6e44f87a2b82.png)
![Earth to Oumouamoua Rendevous](https://user-images.githubusercontent.com/92574647/221439371-51456f0b-f4dc-4cf0-bb50-a7b40c4be592.png)
![Earth to Oumouamoua Fly-by](https://user-images.githubusercontent.com/92574647/221439370-bab4edfd-6ccd-4a95-b552-7a22fb1ce0b8.png)
![Earth to Borisov Fly-by](https://user-images.githubusercontent.com/92574647/221439372-70010d54-6f9c-49f0-af98-20e963eb3b00.png)

## Orbital Element Table
The table below is the orbital elements calculated from the starting radius and velocity vectors for Oumouamoua and Borisov

![Orbital Element Table](https://user-images.githubusercontent.com/92574647/221440324-e1560409-11a3-4edb-8223-ac58d8d3db62.PNG)

Given that the eccentricity of both orbits are greater than one, the trajectory for Oumouamoua and Borisov are interstellar orbits.

I would think that these missions would not be the most realistic to perform due to their margin for error. In most cases, the arrival window for each departure time is very small. If I had to pick a mission, I would pick the rendezvous mission for Borisov. It is enticing to try and pick a different mission with a delta v of 1 or 2 but I think it might be smarter to select a mission with more available arrival dates with acceptable delta V expenditures per each departure dates. This means that the mission could have “hiccups” and still be achievable through secondary thrust inputs to adjust the rendezvous arrival time.
