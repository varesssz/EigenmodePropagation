import warnings
import os.path
import numpy as np
from scipy.io import savemat


class PlaneWaveSet:
    def __init__(
        self,
        k_wave: int | float,  # [rad/m]
        sampling_rate: int | float,  # [1/m]
        window: int | float,  # [m]
    ):
        """
        Calculating plane wave set from the wave vector "k" and the sampling points,
        which is defined by the sampling rate and the (symmetric) window of sampling.
        :param k_wave: [rad/m] Wave vector (k = angular_frequency / propagation_speed = 2 pi / wavelength,
            because k^2 = angular_frequency^2 * permeability * permittivity)
        :param sampling_rate: [1/m] Number of samples taken during one unit (usually given in relative to wavelength).
        :param window: [m] Size of the sampling window.
        """
        self.k_wave = k_wave
        # Numbers of the samples rounded up with int()+1, then making sure it's an odd number
        N = (int(window * sampling_rate) + 1) // 2 * 2 + 1
        # Sampling is made symmetrically for continuous wave likeness
        self.x_sampling = np.linspace(-window/2, window/2, N, endpoint=False)
        # Sampling step size calculated as the mean value of sampling intervals
        self.x_step = np.mean(np.delete(np.delete(np.append(self.x_sampling, 0) - np.append(0, self.x_sampling), -1), 0))
        if self.x_step > (np.pi / k_wave):  # warn if sampling step size is smaller, than half a wavelength
            warnings.warn("Undersampling the expected wavelength!", RuntimeWarning)

        # k_x_max = 2pi/∆x * (N-1)/N as maximum of FFT frequency
        self.k_x_step = (2 * np.pi / self.x_step) * (1 / N)

        # Creating positive k_x[i] array from 0 to k_x_max/2 and negative k_x[i] array from -k_x_max/2 to 0
        self.k_x = np.concatenate((
            np.linspace(
                0,
                self.k_x_step * (N-1) / 2,
                int(N / 2 + 1),  # rounding up half of odd N
            ),
            np.linspace(
                -self.k_x_step * (N-1) / 2,
                0,
                int(N / 2),  # rounding down half of odd N
                endpoint=False,
            )
        ))
        # Throwing away the evanescent plane wave components
        # Also sorting in ascending order
        self.k_x = np.append(
            self.k_x[int(N - k_wave / self.k_x_step + 1):],
            self.k_x[:int(k_wave / self.k_x_step + 1)]
        )
        # Creating k_y from k_x
        self.k_y = np.sqrt(k_wave**2 - self.k_x**2)

        # NumPy handles arctan(np.inf)=pi/2, so zero divide can be ignored
        warnings.filterwarnings("ignore", "divide by zero encountered in divide")
        # Calculating PW incident angles (angle of vector k, measured from the x-axis(CCW))
        self.directions = np.arctan(self.k_y / self.k_x)
        self.directions[self.directions < 0] += np.pi

        # self.param_string = "w%d_x%d_k%d" % (
        #     window,
        #     self.x_sampling.size,
        #     self.k_x.size
        # )

    def save_incident_angles_mat(self, path: str):
        """
        Saving incident angles and sampling points for MATLAB simulation.
        """
        savemat(
            os.path.join(path, "pw_set.mat"),
            dict(
                kwave=self.k_wave,
                sample_points_x=self.x_sampling,
                incident_angles_rad=self.directions,
            ),
        )
