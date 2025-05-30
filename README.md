# Description

This is a Matlab implementation of Wu et al. that is being developed for the following reasons:

- No public implementation of this model exists (to our knowledge)
- We would like to be able to reproduce the plots published in the paper

Wu JZ, Herzog W, Cole GK. Modeling dynamic contraction of muscle using the cross-bridge theory. Mathematical biosciences. 1997 Jan 1;139(1):69-78. https://doi.org/10.1016/S0025-5564(96)00115-0


# Quick start guide:

1. Clone this repository
2. Start Matlab
3. Run main_HuxleyModel.m

# Status

The solutions produced have oscillations in the spatial domain. The implementation of the method of lines has been tested against the PDE that describes beam bending, and this is fine. At the moment we are digging into the details to see if we can isolate this problem as no such spatial oscillation appears in Wu et al.




