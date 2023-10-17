import numpy as np
import pickle

from MyPlotlyFigure import MyPlotlyFigure
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

    print("Load Transfer-matrix from file...")
    with open("data/transfer_mat_w120_x5605_k801_structA.pkl", "rb") as file:
        transfer_matrix = pickle.load(file)

    # Create excitation vector
    excitation_vector = np.zeros(pw_set.k_x.size, np.complex64)
    excitation_vector[469] = 1
    excitation_vector[400] = 1
    excitation_vector[93] = 1

    # Calculate Ez(k_x) with transfer-matrix and excitation vector
    print("Calculate resulted plane waves from transfer-matrix...")
    ez_kx = transfer_matrix @ excitation_vector

    # Reconstruct total E-field with this result
    y_reconstruct = 0.0
    print("Reconstruct total E-field from plane wave components at y=%.1fm..." % y_reconstruct)
    ez_reconstructed_1d = np.zeros(pw_set.x_sampling.size, dtype=np.complex64)
    for k_x_i, k_y_i, ez_kx_forward_i in zip(pw_set.k_x, pw_set.k_y, ez_kx):
        ez_reconstructed_1d += (
            ez_kx_forward_i * np.exp(1j * k_y_i * y_reconstruct) * np.exp(1j * k_x_i * pw_set.x_sampling)
        )

    # Load in simulation results for the same excitation
    print("Load previously simulated validation data from file...")
    e_z_simulated = np.genfromtxt(
        fname="data/e_field_y0_w120_x5605_k801_structA.txt",
        dtype=np.complex64,
        delimiter=",",
    )

    # Adding reconstruction and validation to total E-field plot for comparison
    print("Create Plotly figure with the reconstructed and the validation total E-field plots...")
    fig = MyPlotlyFigure()
    fig.update_layout(
        xaxis_title="x [m]",
        yaxis_title="E_z real value [V/m]",
    )
    fig.add_scatter(
        x=pw_set.x_sampling,
        y=np.real(ez_reconstructed_1d),
        name=f"reconstructed (y = {y_reconstruct})",
    )
    fig.add_scatter(
        x=pw_set.x_sampling,
        y=np.real(e_z_simulated),
        name=f"simulated (y = {y_reconstruct})",
    )
    fig.matlab_styling()
    fig.show()

    # Root-Mean-Square Deviation
    print("Calculate RMS Deviation of simulated and reconstructed E-field real values...")
    rmsd = (np.real(e_z_simulated) - np.real(ez_reconstructed_1d))**2
    rmsd = np.sqrt(rmsd.sum(axis=0) / rmsd.size)
    nrmsd = rmsd / (np.real(e_z_simulated).max(axis=0) - np.real(e_z_simulated).min(axis=0))
    print(f"\nRMSD:{np.format_float_scientific(rmsd, 2)},   NRMSD:{np.format_float_scientific(nrmsd, 2)}")
