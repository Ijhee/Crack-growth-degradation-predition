<div align=center> <h1> Crack-growth-degradation-predition  </h1> </div>
<div align=center> <h4> 24-1 Prognostics and Health Management  </h4> </div>

## 1. Problem Defintion
<p align="center"><img src=https://github.com/Ijhee/Crack-growth-degradation-predition/assets/96717686/64f36e3d-3edd-464c-af45-497d9c533f80></p>
- Using three methods—nonlinear regression, Bayesian models, and particle filters to predict crack growth</br>
- The degradation data used for the analysis is shown in the figure above. </br>
<p align="center"><img src=https://github.com/Ijhee/Crack-growth-degradation-predition/assets/96717686/827653e7-ef0d-42d3-90f0-66656e5c7d10S></p>
- Fig. 2 shows the code for problem definition in nonlinear regression, which includes defining the data to be analyzed, the threshold, and the true values.
<p align="center"><img src=https://github.com/Ijhee/Crack-growth-degradation-predition/assets/96717686/714144fa-d23e-4cec-b335-26adc3121d29></p>
- Fig. 3 shows the problem definition code for the Bayesian model, where the problem definition is the same as in Fig. 2, but with the addition of the 'initDisPar' matrix for the parameters of the initial/prior distribution and the MCMC-specific parameter 'burnIn'.
<p align="center"><img src=https://github.com/Ijhee/Crack-growth-degradation-predition/assets/96717686/4cafb018-0a67-4d32-bad9-898449625063></p>
- Fig. 4 shows the problem definition code for the particle filter, where the problem definition is the same as in Fig. 3, but since the initial crack is assumed to be 0.01m, the standard deviation is also set to 0.
