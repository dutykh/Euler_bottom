# A full Euler solver on general bottoms

The present Matlab code is an implementation of the full Euler equations solver based on the method of conformal variables. The pecularity here is that the solver works on general (but smooth) bottoms. The method is described in the reference given below. In a few words it is a Fourier-type pseudo-spectral solver. Standard Matlab time stepper is used to advance the solution in time. The solution is expected to be spectrally accurate. See the code to set the appropriate error tolerance parameters.

Any comments/suggestions/bug reports are welcome!

To contact me, please, consult my web page:

* [www.denys-dutykh.com](http://www.denys-dutykh.com/)

![Snapshot](pics/InitialSnapshot.png)

## Reference:

* C. Viotti, **D. Dutykh**, F. Dias. *The conformal-mapping method for surface gravity waves in the presence of variable bathymetry and mean current*, Procedia IUTAM, **11**, 110 - 118, 2014

## Acknowledgements:

The author would like to thank Professor [Yury Stepanyants](https://staffprofile.usq.edu.au/profile/yury-stepanyants) (University of Southern Queensland, Australia) for bringing my attention to this problem.

---
