#!/bin/bash

# Author: Radu Cimpeanu
# Date: 11/08/2024

# Additional velocities or resolution levels can be added below
for V0 in 0.9392; do
	for LEVEL in 12 13 14; do

		# Copy all files to renamed folder based on key parameters
		cp -r MasterImpactSingleVOF/ Impact-V$V0-Level$LEVEL-SingleVOF
		cd Impact-V$V0-Level$LEVEL-SingleVOF/
		
		# Compile code to create the executable (including visualisation)
		qcc -O2 -w -fopenmp -Wall DropImpact.c -lm -o DropImpact -L$BASILISK/gl -lglutils -lfb_tiny -lGLU -lGLEW -lGL -lX11

		# Specify parallelisation features		
		export OMP_NUM_THREADS=4
		
		# parameters:
		# 1. rho liquid (dimensional)
		# 2. rho gas (dimensional)
		# 3. mu liquid (dimensional)
		# 4. mu gas (dimensional)
		# 5. sigma (surface tension coefficient, dimensional)
		# 6. g acceleration (dimensional)
		# 7. drop radius (dimensional)
		# 8. initial drop velocity (dimensional)
		# 9. simulation end time (dimensionless, tailored to contact and bounce duration)
		# 10. max level

		# Run executable
		./DropImpact 873.0 1.21 1.7e-3 1.81e-5 0.0187 9.81 0.17e-3 $V0 3.001 $LEVEL
		
		cd ..
	done
done 
