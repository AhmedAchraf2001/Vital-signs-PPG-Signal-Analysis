# PPG Signal Peak detection using python 

**Description :**  Systolic upslopes are detected from a signal generated with a slope sum function, which sums the magnitudes of the PPG upslopes in the previous 0.17 s. Adaptive thresholding is used to identify systolic upslopes in this signal. The 'qppgfast' implementation of this beat detector was used, after testing showed it performed similarly to the original 'qppg' implementation.


**inputs**

	    sig : a vector of PPG values

	    fs : the sampling frequency of the PPG in Hz

**outputs**

	    peaks : indices of detected pulse peaks

	    onsets : indices of detected pulse onsets

## My Contribuation
my contribution was about implementation the algorithm in Python programming language.


## Reference

A. N. Vest et al., 'An open source benchmarked toolbox for cardiovascular waveform and interval analysis,' Physiological Measurement, vol. 39, no. 10, 2018. [https://doi.org/10.1088/1361-6579/aae021](https://doi.org/10.1088/1361-6579/aae021)


## Author

-   Several authors have contributed to the code (see below)
    
-   Peter H. Charlton - did very little, just wrote this wrapper (which detects peaks given the onsets provided by the code)
