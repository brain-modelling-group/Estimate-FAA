# ECG-to-FAA


**Scripts for ECG to FAA analysis**
 - training_MDL_hrv.mat: contains the gaussian process regression model from "Bedside tracking of functional autonomic age in preterm infants" manuscript 
 - ecg_to_faa_example: a brief run through of ECG epoch selection (example ECG epoch provided; option to reject/accept epochs), calculation of HRV features from the ECG and age estimation based on trained model 
 - ecg_to_nn_estimation: code to extract the NN interval from the ECG
 - calculate_features: the HRV features to be extracted following NN interval estimation

**Dependencies**

MATLAB, tested on versions 2020b and 2022a


**Citation**

If our scripts are used to inform your work, please cite us: 

Iyer, K.K., Leitner, U., Giordano, V. et al. Bedside tracking of functional autonomic age in preterm infants. Pediatr Res (2022). https://doi.org/10.1038/s41390-022-02376-2
