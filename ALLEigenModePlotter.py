import numpy as np
import pickle
from MyPlotlyFigure import MyPlotlyFigure


if __name__ == '__main__':

    # pw_set_model = "plane_wave"
    # pw_set_model = "points_along_line"
    # pw_set_model = "points_along_circle"
    pw_set_model = "points_as_phased_array"

    e_field_all_eigen = np.genfromtxt(
        fname="MATLAB/output/all_eigenmodes_by_%s.txt" % pw_set_model,
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
    )
    fig.update_layout(
        xaxis_title="x [m]",
        # yaxis_title="eigenmode indexed by descending eigenvalue order",
        yaxis_range=[0, 800]
    )
    fig.show()

    matlab_structure = "structureC"

    print("Load Transfer-matrix from file...")
    with open("data/transfer_mat_%s_from_%s.pkl" % (matlab_structure, pw_set_model), "rb") as file:
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
        yaxis_range=[0, 800]
    )
    fig.show()
