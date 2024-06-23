<div align=center> <h1> Crack-growth-degradation-predition  </h1> </div>
<div align=center> <h4> 24-1 Prognostics and Health Management  </h4> </div>

## 1. Problem Defintion
### 1.1 Problem Definition
<p align="center"><img src=https://github.com/Ijhee/Crack-growth-degradation-predition/assets/96717686/64f36e3d-3edd-464c-af45-497d9c533f80></p>

- Using three methods—nonlinear regression, Bayesian models, and particle filters to predict crack growth.
- The degradation data used for the analysis is shown in the figure above. </br>
<p align="center"><img src=https://github.com/Ijhee/Crack-growth-degradation-predition/assets/96717686/827653e7-ef0d-42d3-90f0-66656e5c7d10S></p>

- Fig. 2 shows the code for problem definition in nonlinear regression, which includes defining the data to be analyzed, the threshold, and the true values.
  
<p align="center"><img src=https://github.com/Ijhee/Crack-growth-degradation-predition/assets/96717686/714144fa-d23e-4cec-b335-26adc3121d29></p>

- Fig. 3 shows the problem definition code for the Bayesian model, where the problem definition is the same as in Fig. 2, but with the addition of the 'initDisPar' matrix for the parameters of the initial/prior distribution and the MCMC-specific parameter 'burnIn'.
  
<p align="center"><img src=https://github.com/Ijhee/Crack-growth-degradation-predition/assets/96717686/4cafb018-0a67-4d32-bad9-898449625063></p>

- Fig. 4 shows the problem definition code for the particle filter, where the problem definition is the same as in Fig. 3, but since the initial crack is assumed to be 0.01m, the standard deviation is also set to 0.
  
### 1.2 Degradation Model

<p align="center"><img src=https://github.com/Ijhee/Crack-growth-degradation-predition/assets/96717686/b39637cb-9d41-49d5-8ef3-7c3be9ed0a63></p>

- Nonlinear least square and Bayesian models used Eq. 1 as the degradation model, while the Particle filter used Eq. 2.
  
## 2. Modifying the Codes for the Crack Growth Example & Results

<p align="center"><img src=https://github.com/Ijhee/Crack-growth-degradation-predition/assets/96717686/8c4a359a-bae2-4988-915a-b02bee76aa18></p>
<p align="center"><img src=https://github.com/Ijhee/Crack-growth-degradation-predition/assets/96717686/aaed4523-5b16-400f-8d74-5705f0f7f01c></p>
<p align="center"><img src=https://github.com/Ijhee/Crack-growth-degradation-predition/assets/96717686/4b421bcc-ec53-4ec7-819f-cc29960b237e></p>

- The degradation prediction results using nonlinear regression, Bayesian models, and particle filters are shown in Figures 5, 6, and 7.
- From these figures, it is evident that the results of NLS (nonlinear least squares) exhibit greater uncertainty compared to the other methods.
- This is because nonlinear regression cannot utilize prior information, leading to parameter estimates over a wide range, which can cause significant variations in degradation predictions.
- In contrast, Bayesian models and particle filters use prior information, resulting in parameters being identified within a narrower range compared to nonlinear regression.
- Consequently, the prediction uncertainty is relatively less than that of nonlinear regression.
