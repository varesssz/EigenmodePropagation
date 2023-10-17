import warnings
import numpy as np
from scipy.io import savemat


def create_pw_set(
        k_wave: int | float,  # [rad/m]
        sampling_rate: int | float,  # [1/m]
        window: int | float,  # [m]
        save_incident_angles_mat=False,
):
    """
    Calculating plane wave set from the wave vector "k" and the sampling points, which is defined by the sampling rate
    and the (symmetric) window of sampling. The incident wave angles also calculated and saved for MATLAB usage.
    :param k_wave: [rad/m] Wave vector (k = angular_frequency / propagation_speed = 2 pi / wavelength,
        because k^2 = angular_frequency^2 * permeability * permittivity)
    :param sampling_rate: [1/m] Number of samples taken during one unit (usually given in relative to the wavelength).
    :param window: [m] Size of the sampling window.
    :param save_incident_angles_mat: [bool] To save incident angles (and sampling points) to MATLAB file (.mat) or not.
    :return: Dictionary, containing the sampling points ["x_sampling"], the sampling step size ["x_step"],
        the x and y components of the wave vectors ["k_x"] and ["k_y"], and the step size of k_x ["k_x_step"].
    """
    # Numbers of the samples rounded up with int()+1, then making sure it's an odd number
    N = (int(window * sampling_rate) + 1) // 2 * 2 + 1
    # Sampling is made symmetrically for continuous wave likeness
    x_sampling = np.linspace(-window/2, window/2, N, endpoint=False)
    # Sampling step size calculated as the mean value of sampling intervals
    x_step = np.mean(np.delete(np.delete(np.append(x_sampling, 0) - np.append(0, x_sampling), -1), 0))
    if x_step > (np.pi / k_wave):  # warn if sampling step size is smaller, than half a wavelength
        warnings.warn("Undersampling the expected wavelength!", RuntimeWarning)

    # k_x_max = 2pi/∆x * (N-1)/N as maximum of FFT frequency
    k_x_step = (2 * np.pi / x_step) * (1 / N)

    # Creating positive k_x[i] array from 0 to k_x_max/2 and negative k_x[i] array from -k_x_max/2 to 0
    k_x = np.concatenate((
        np.linspace(
            0,
            k_x_step * (N-1) / 2,
            int(N / 2 + 1),  # rounding up half of odd N
        ),
        np.linspace(
            -k_x_step * (N-1) / 2,
            0,
            int(N / 2),  # rounding down half of odd N
            endpoint=False,
        )
    ))
    # Throwing away the evanescent plane wave components
    # Also sorting in ascending order
    k_x = np.append(
        k_x[int(N - k_wave / k_x_step + 1):],
        k_x[:int(k_wave / k_x_step + 1)]
    )
    # Creating k_y from k_x
    k_y = np.sqrt(k_wave**2 - k_x**2)

    # Calculating PW incident angles (angle of vector k, measured from the x-axis(CCW))
    warnings.filterwarnings("ignore", "divide by zero encountered in divide")  # arctan(np.inf)=pi/2 can be handled
    incident_angles = np.arctan(k_y / k_x)
    incident_angles[incident_angles < 0] += np.pi

    # Saving incident angles and sampling points for MATLAB simulation
    if save_incident_angles_mat:
        savemat(
            "./data/pw_set_w%d_x%d_k%d.mat" % (window, x_sampling.size, incident_angles.size),
            dict(
                kwave=k_wave,
                sample_points_x=x_sampling,
                incident_angles_rad=incident_angles,
            ),
        )

    return dict(
        x_sampling=x_sampling,
        x_step=x_step,
        k_x=k_x,
        k_x_step=k_x_step,
        k_y=k_y,
    )
