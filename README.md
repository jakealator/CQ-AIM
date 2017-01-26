# Wave equation simulation with convolution quadrature and adaptive integral method

Simulates propogation of acoustic waves through an inhomogeneous medium with 
shape and material properties specified by the user. Note that code is currently 
still very much in progress.

Specifically, solves the time domain Lippmann-Schwinger equation using a 
convolution quadrature time stepper and sped-up piecewise constant Galerkin 
spatial discretization. Speed up is via the adaptive integral methid (AIM), 
which computes near field contributions with a Galerkin approach and far 
field contributions with a fast Fourier transform (FFT) approach. For many
more details on this, see chapter 5 of my Ph.D thesis (which is not yet publically available).  

## Getting Started

To run this code, download all files in the repository to the same folder. 
Run Matlab through that folder.

