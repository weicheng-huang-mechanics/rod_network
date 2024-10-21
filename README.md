# Overview

This project focused on the 3d rod network dynamics.

<br/><img src='demo2.gif' width="600">

To run this code, you should have a Linux Ubuntu system

# Make
g++ -I /usr/local/include/eigen-3.3.7/ main.cpp world.cpp setInput.cpp timeStepper.cpp inertialForce.cpp externalGravityForce.cpp dampingForce.cpp elasticStretchingForce.cpp elasticBendingForce.cpp elasticTwistingForce.cpp elasticPlate.cpp -lGL -lglut -lGLU -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -llapack -lgfortran -fopenmp -lpthread -lm -Ofast -o simDER

# Run
./simDER option.txt
