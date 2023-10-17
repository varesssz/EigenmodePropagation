import numpy as np
import pickle

import AskUserinputRecursively
from PlaneWaveSet import PlaneWaveSet

if __name__ == '__main__':

    # Parameters of PWD
    print("Initialization...")
    wavelength = 299792458 / 1e9  # [m]
    k_wave = 2 * np.pi / wavelength
    # Initialize plane wave set
    pw_set = PlaneWaveSet(
        k_wave=k_wave,
        sampling_rate=14 / wavelength,
        window=120,
    )
    pw_set.save_incident_angles_mat()

    # Only go forward if MATLAB simulation is done and the E-field transfer-matrix is ready
    if not AskUserinputRecursively.yes_or_no("Continue with the process of the MATLAB simulation results?"):
        raise KeyboardInterrupt
    print("Reading MATLAB results...")
    matlab_structure = "structA"
    simulation_string = "w%d_x%d_k%d" % (
        np.round(pw_set.x_sampling.size * pw_set.x_step),
        pw_set.x_sampling.size,
        pw_set.k_x.size
    )
    efield_transfer_matrix = np.genfromtxt(
        fname="data/efield_transfer_mat_%s_%s.txt" % (simulation_string, matlab_structure),
        dtype=np.complex64,
        delimiter=",",
    )
    cylinder_pos = np.genfromtxt(
        fname="data/cylinders_%s.txt" % matlab_structure,
        dtype=np.float32,
        delimiter=",",
    )

    # Plane Wave Decomposition with DFT
    print("Plane Wave Decomposition with DFT...")
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
    ez_kx = np.transpose(ez_kx)

    # Save resulting Transfer-matrix in file
    print("Saving transfer-matrix...")
    with open("data/transfer_mat_%s_%s.pkl" % (simulation_string, matlab_structure), "wb") as file:
        pickle.dump(ez_kx, file)

    print("DONE")
