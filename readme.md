Kuramoto-Sakaguchi Model Simulation Applet by Josu Blanco and Elena del Campo. 2022.

This applet simulates the evolution of a group of N coupled oscillators according to the Kuramoto-Sakaguchi model.

## Required libraries
- Python 3
- Numpy, matplotlib and PyQt5

## Installation
The downloaded folder consists of 4 .py Python3 and 2 .txt files. To execute the applet just run KSModel.py

Windows:
1. Open the terminal by typing cmd on the searchbar.
2. Go to the directory where KSModel.py is by typing cd C:\...\KSModel.py (change ... for the correct path).
3. Type python KSModel.py

Linux:
1. Open the terminal
2. Go to the directory where KSModel.py is by typing cd /.../KSModel.py (change ... for the correct path).
3. Run the file with ./KSModel.py

## Usage
The applet is divided into 4 tabs. 

On the far left we find the current state of the system. Each point represents the angular motion of an individual oscillator. The colour of the dot represents the intrinsic angular frequency of that oscillator. There is also a small "worm" alongside a cross in the centre. The distance from the worm to the centre represents the modulus of the order parameter in the current state of the system. The position of the worm represents the value of the mean angle. The worm has a length of 20 time steps so that the change of the system in the recent past can be seen.

On the top-center we find the Gaussian frequency distribution graph. This graph consists of a histogram that scales with N and the proper gaussian function for the given Deviation. The colours of the histogram (discrete) relate to the colours of the dots in the state of the system window (continous colour map).

On the bottom-center we find the modulus of the order parameter vs. the elapsed time steps. This live graph is joined by a static dotted red line corresponding to the mean of the fluctuations of the system, for coupling constant equals zero, given by N.

On the right side we find the interactive control boxes, these boxes are coded with a two-colour system. 
Boxes in colour light-red can be changed live and correspond to:
- Slow time: Slows the graphing time, where 0 is fastest and 100 slowest.
- Number of osc. on screen: Random slice of the total number of oscillators for graphing.
- Coupling constant: Change the coupling constant K.
- Worm: Tickbox for deleting the worm in the current system tab.

Boxes in colour light-blue are changed only on reset and correspond to:
- Number of oscillators: Total number of oscillators for calculation.
- Iteration period: Changes the calculation period. A period of 2 means that the program will calculate the ODE every 2 time units.
- Deviation: Changes the deviation sigma of the gaussian distribution. This is linked to the critical value for K.
- Alpha: Changes the phasing inside the sine.
- Close: Closes the program.
- Reset: Resets the program saving the current state of the control boxes to the new_vars.txt file.

