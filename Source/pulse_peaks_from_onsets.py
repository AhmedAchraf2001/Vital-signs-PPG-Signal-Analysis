def pulse_peaks_from_onsets(sig, onsets):
    import numpy as np 
    peaks = []
    for wave_no in range(1, len(onsets)-1):
        temp = np.argmax( sig[onsets[ wave_no ] : onsets[ wave_no + 1] ] )
        peaks.append(temp + onsets[wave_no] - 1)
    return peaks