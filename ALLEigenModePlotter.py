import numpy as np
import pickle
from MyPlotlyFigure import MyPlotlyFigure
from PlaneWaveSet import PlaneWaveSet


if __name__ == '__main__':

    e_field_all_eigen = np.genfromtxt(
        fname="MATLAB/data/structC_w120_x5605_k801/0deg_to_180deg_excitation/eigenmode_y=-1_5_all.txt",
        dtype=np.complex64,
        delimiter=",",
    )

    fig = MyPlotlyFigure()
    fig.matlab_styling()
    fig.add_heatmap(
        x=np.linspace(-4.5, 4.5, 60 + 1),
        z=np.abs(e_field_all_eigen),
        colorbar_title_text="E_z absolute value [V/m]",
        colorbar_title_side="right",
        # colorscale="jet",
        # zmid=0,
        # zmax=10,
        # zmin=-10,
    )
    fig.update_layout(
        xaxis_title="x [m]",
        # yaxis_title="eigenmode indexed by descending eigenvalue order",
        yaxis_range=[720, 800]
    )
    fig.show()

    print("Initializing plane wave set...")
    wavelength = 299792458 / 1e9  # [m]
    k_wave = 2 * np.pi / wavelength
    pw_set = PlaneWaveSet(
        k_wave=k_wave,
        sampling_rate=14 / wavelength,
        window=120,
    )

    matlab_structure = "structureC"

    print("Load Transfer-matrix from file...")
    with open("data/transfer_mat_%s.pkl" % matlab_structure, "rb") as file:
        transfer_matrix = pickle.load(file)

    # Eigen values and vectors of the transfer matrix
    eig_values, eig_vectors = np.linalg.eig(transfer_matrix)

    # Sort eigen values and vectors by eigen value
    sort = np.flip(np.argsort(np.abs(eig_values)))
    eig_values = eig_values[sort]
    eig_vectors = eig_vectors[:, sort]

    fig.data = []
    fig.add_scatter(
        x=np.abs(eig_values),
        mode="markers",
        marker_size=10,
    )
    fig.update_layout(
        xaxis_title="$|~\lambda~|$",
        xaxis_range=[0, 1.05],
        yaxis_title="eigenmode indices (descending eigenvalue order)",
        yaxis_range=[720, 800]
    )
    fig.show()
