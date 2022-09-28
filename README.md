# Bouncing droplets
Direct numerical simulation code infrastructure for drop impact onto liquid pools in the inertio-capillary regime, supporting collaborative work with the Harris Lab at Brown.  

Complements the preprint available at https://arxiv.org/abs/2209.13276 mathematical modelling repository at https://github.com/harrislab-brown/BouncingDroplets.

## Installation
* The code relies on [Basilisk](<http://basilisk.fr/>) to model the Navier-Stokes equations. See the [installation page](<http://basilisk.fr/src/INSTALL>) for instructions. 
* Full visualisation capabilities have been used in order to generate animations. These may be switched off depending on the local architecture.
* The two-phase non-coalescing fluid volume implementation by V. Sanjay available [here](https://github.com/VatsalSy/Lifting-a-sessile-drop/blob/master/CaseI/two-phaseDOD.h) has been successfully employed in this study to limit numerical artifacts during contact time.

## Running the code
Once the Basilisk structure is in place, the driver code here is built in order to navigate parameter sweeps in velocity $V_0$ and resolution level, with one of each values added to the run_master_example.sh for brevity. Other parameters can be varied through this shell script, with both physical and computational handles provided. 

The code can be executed by simply executing this shell script via *sh run_master_example.sh* inside a terminal. Output will then be produced within a foldering structure that consists of summary DNS execution information, mass conservation and VOF data, interface coordinates, simulation slices and animations, which can be used for further post-processing.

## Example results
The uploaded framework provides a subset of the data generated for the test case described as part of Figure 3 in [the accompanying manuscript](https://arxiv.org/abs/2209.13276) and represents the case of a water droplet impinging onto a liquid pool at a moderate impact velocity $V_0 = 0.3855$ m/s.
