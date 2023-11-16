import warnings
import os.path
import numpy as np
from scipy.io import savemat
from scipy.signal import windows


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
        self.x_step = np.mean(np.diff(self.x_sampling))
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

        self.model = "plane_wave"

        self.parameters_to_save = dict(
            kwave=self.k_wave,
            sample_points_x=self.x_sampling,
            inc_angles=self.directions,
            model_select=self.model,
        )

    def set_up_points_along_line_model(
            self,
            array_distance: int | float,
    ):
        self.model = "points_along_line"
        # Source position from line position and direction
        sources_y = array_distance
        sources_x = sources_y / np.tan(self.directions)
        sources_xy = sources_x + 1j * sources_y

        # Save it to parameter dictionary
        self.parameters_to_save["sources_pos"] = sources_xy
        self.parameters_to_save["model_select"] = self.model

    def set_up_points_along_circle_model(
            self,
            array_radius: int | float,
    ):
        self.model = "points_along_circle"
        # Source position from line position and direction
        sources_x = np.cos(self.directions) * array_radius
        sources_y = np.sin(self.directions) * array_radius
        sources_xy = sources_x + 1j * sources_y

        # Save it to parameter dictionary
        self.parameters_to_save["sources_pos"] = sources_xy
        self.parameters_to_save["model_select"] = self.model

    def set_up_phased_array_model(
            self,
            array_distance: int | float,
            array_length: int | float,
            element_dist_per_lambda: int | float,
            taylor_windowing=False,
    ):
        self.model = "points_as_phased_array"
        # Parameters of the model
        element_distance = element_dist_per_lambda * (2 * np.pi / self.k_wave)
        element_count = (int(array_length / element_distance) + 1) // 2 * 2 + 1

        # Source position from line position and direction
        sources_x = np.linspace(-array_length / 2, array_length / 2, element_count)
        sources_y = array_distance
        sources_xy = sources_x + 1j * sources_y

        # Weighting antenna elements
        weights = np.ones(element_count)
        if taylor_windowing:
            weights = windows.taylor(element_count)

        # Save it to parameter dictionary
        self.parameters_to_save["sources_pos"] = sources_xy
        self.parameters_to_save["sources_w"] = weights
        self.parameters_to_save["model_select"] = self.model

    def save_parameters_for_matlab(self, path: str):
        """
        Saving incident angles, sampling points and model parameters for MATLAB simulation.
        """
        savemat(
            file_name=os.path.join(path, "pw_set.mat"),
            mdict=self.parameters_to_save,
        )
