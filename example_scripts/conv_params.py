
# aw 2025-01-25 12:16:01 - Updated to include mega_params
conv_params = {
    'binSize': 0.01,
    'gaussianSigma': 0.01,
    'thresholdBurst': 1.0,
    'min_peak_distance': None, # no minimum peak distance - given the prominence method used in the burst detection, this should be fine...I hope
    'prominence': 2,
}

# TODO: this isn't implemented everywhere yet... verify. # aw 2025-01-25 12:18:58 - but things should work without it for now. Hardcoded into some functions.
mega_params = { 
    'binSize': conv_params['binSize']*3,
    #'gaussianSigma': conv_params['gaussianSigma']*15,
    ## aw 2025-02-04 17:33:19 - I want to be a little more sensitive to the peaks in the mega data
    'gaussianSigma': conv_params['gaussianSigma']*8,
    'thresholdBurst': 1.0,
    'min_peak_distance': None, 
    'prominence': 2,
}
