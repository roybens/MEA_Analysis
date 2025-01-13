import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline



def polynomial_interpolation(trace, start_idx, end_idx, degree=3):
    # Define x and y points outside the artifact region
    x = np.concatenate((
        np.arange(0, start_idx),  # Points before
        np.arange(end_idx + 1, len(trace) - 1)   # Points after
    ))
    y = np.concatenate((
        trace[0: start_idx],
        trace[end_idx + 1 : len(trace) - 1]
    ))

    # Fit polynomial to surrounding points
    # poly_coeffs = np.polyfit(x, y, degree)

    
    spline = CubicSpline(x, y)


    # Interpolate the artifact region
    artifact_indices = np.arange(start_idx, end_idx + 1)
    # trace[artifact_indices] = np.polyval(poly_coeffs, artifact_indices)
    trace[artifact_indices] = spline(artifact_indices)

    # plt.plot(x, y, 'o', label="Surrounding Points")
    # plt.plot(artifact_indices, trace[artifact_indices], 'r-', label="Interpolated Region")
    # plt.legend()
    # plt.show()

    return trace