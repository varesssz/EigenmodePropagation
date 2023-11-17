import numpy as np
import pickle
from scipy.io import savemat

import AskUserinputRecursively
from MatlabRunner import MatlabRunner
from PlaneWaveSet import PlaneWaveSet

if __name__ == '__main__':

    print(">> Initializing plane wave set...")
    wavelength = 299792458 / 1e9  # [m]
    k_wave = 2 * np.pi / wavelength
    # Initialize plane wave set
    pw_set = PlaneWaveSet(
        k_wave=k_wave,
        sampling_rate=14 / wavelength,
        window=120,
    )

    # Set up the simulation model
    # Plane wave model is used if no setup is uncommented from below

    # pw_set.set_up_points_along_line_model(
    #     array_distance=-2.0 - 10.0,
    # )

    # pw_set.set_up_points_along_circle_model(
    #     array_radius=-2.0 - 10.0,
    # )

    pw_set.set_up_phased_array_model(
        array_distance=-2.0 - 3.0,
        array_length=20,
        element_dist_per_lambda=0.5,
        taylor_windowing=False,
    )

    # Save PW set's parameters for MATLAB usage
    pw_set.save_parameters_for_matlab(path="./MATLAB/data/")

    cylinder_structure_name = "structureC"
    savemat("./MATLAB/data/cylinder_struct.mat", dict(
        cylinders=np.genfromtxt(
            fname="data/cylinders_%s.txt" % cylinder_structure_name,
            dtype=np.float64,
            delimiter=",",
            skip_header=2,
        )
    ))

    if AskUserinputRecursively.yes_or_no("Run MATLAB simulation to excite with every PW in the set?"):
        print(">> Simulating with the model %s" % pw_set.model)
        matlab = MatlabRunner()
        matlab.run_matlab_script("mieScatt_multiPW.m")
        matlab.export_fixer()

    print(">> Reading MATLAB results...")
    efield_transfer_matrix = np.genfromtxt(
        fname="MATLAB/data/efield_transfer_mat.txt",
        dtype=np.complex64,
        delimiter=",",
    )

    # Plane Wave Decomposition with DFT
    print(">> Performing Plane Wave Decomposition with DFT...")
    # Allocate memory for plane waves (E_z(k_x))
    ez_kx = np.empty((pw_set.k_x.size, pw_set.k_x.size), np.complex64)
    # Calculate PWD for each column (for loop goes through rows, so transposing is necessary)
    for i, e_field in enumerate(np.transpose(efield_transfer_matrix)):
        # E_z(k_x, y) with FFT{E_z(x, y)}
        ez_kx_y = np.fft.fft(e_field) / pw_set.x_sampling.size

        # Throwing away results related to the evanescent plane wave components
        # Also sorting in ascending order
        ez_kx_y = np.append(
            ez_kx_y[int(pw_set.x_sampling.size - k_wave / pw_set.k_x_step + 1):],
            ez_kx_y[:int(k_wave / pw_set.k_x_step + 1)]
        )

        # Compensating that the sampling is not starting from x=0 coordinate
        # Saving it in the E_z(k_x) array
        ez_kx[i] = ez_kx_y * np.exp(1j * pw_set.k_x * pw_set.x_sampling[0])

        # Compensation with e^(j k_y y) that the sampling is not done through the origin is not necessary
        # As the simulation sampling was done at the right place (y=0) where the phase center of the PWs are
        # ez_kx[i] = ez_kx[i] * np.exp(1j * k_y * 0)

    # Transpose back, so that the k_x variable goes along the columns
    transfer_matrix = np.transpose(ez_kx)

    # Save resulting Transfer-matrix in file
    print(">> Saving transfer-matrix...")
    with open("data/transfer_mat_%s_from_%s.pkl" % (cylinder_structure_name, pw_set.model), "wb") as file:
        pickle.dump(transfer_matrix, file)
