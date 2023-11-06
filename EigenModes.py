import numpy as np
import pickle
from scipy.io import savemat

import AskUserinputRecursively
from MatlabRunner import MatlabRunner
from PlaneWaveSet import PlaneWaveSet
from MyPlotlyFigure import MyPlotlyFigure


def eigenmode_filter(
        eigen_values,
        eigen_vectors,
        pw_directions,
        filter_eig_abs_value_in=(0.0, 1.0),
        filter_eig_angle_in=(-180.0, 180.0),
        filter_eig_vector_above=0.00,
        filter_eig_vector_only_in_direction=(0, 180),
        filter_eig_vector_direction_hitrate=0.0,
):
    # Filter eigen values that falls in between of "filter_eig_value_in"
    boolarr_eig_value = np.bitwise_and(
        np.bitwise_and(
            filter_eig_abs_value_in[0] <= np.abs(eigen_values),
            np.abs(eigen_values) <= filter_eig_abs_value_in[1]
        ),
        np.bitwise_and(
            filter_eig_angle_in[0] <= np.rad2deg(np.angle(eigen_values)),
            np.rad2deg(np.angle(eigen_values)) <= filter_eig_angle_in[1]
        )
    )
    # Investigate eigen vector PW components in given incident angle range "filter_eig_vector_in_incident_angles"
    boolarr_directions = np.bitwise_and(
        filter_eig_vector_only_in_direction[0] <= np.rad2deg(pw_directions),
        np.rad2deg(pw_directions) <= filter_eig_vector_only_in_direction[1]
    )
    # Filter eigen vectors that is above given value in the given incident angle range
    boolmat_eig_vector = filter_eig_vector_above < np.abs(eigen_vectors[boolarr_directions.nonzero()])
    # Filter eigen vectors if previous statement is True for at least XX% of those angles
    boolarr_eig_vector = filter_eig_vector_direction_hitrate * boolarr_directions.sum() <= boolmat_eig_vector.sum(0)

    return np.bitwise_and(boolarr_eig_value, boolarr_eig_vector).nonzero()[0]


if __name__ == '__main__':

    print("Initializing plane wave set...")
    wavelength = 299792458 / 1e9  # [m]
    k_wave = 2 * np.pi / wavelength
    pw_set = PlaneWaveSet(
        k_wave=k_wave,
        sampling_rate=14 / wavelength,
        window=120,
    )

    matlab_structure = "structureB"

    print("Load Transfer-matrix from file...")
    with open("data/transfer_mat_%s.pkl" % matlab_structure, "rb") as file:
        transfer_matrix = pickle.load(file)

    # Eigen values and vectors of the transfer matrix
    eig_values, eig_vectors = np.linalg.eig(transfer_matrix)

    # Sort eigen values and vectors by eigen value
    sort = np.flip(np.argsort(np.abs(eig_values)))
    eig_values = eig_values[sort]
    eig_vectors = eig_vectors[:, sort]
    # Save all eigen vectors for "mieScatt_ALLeigenPW.m" script to simulate with all
    savemat(
        file_name="MATLAB/data/eigen_vectors.mat",
        mdict=dict(
            eigen_vectors=eig_vectors,
        )
    )

    # Select eigenmodes from the eigenvectors by filtering eigenvalues and eigenvector elements
    indices = eigenmode_filter(
        eig_values,
        eig_vectors,
        pw_set.directions,
        filter_eig_abs_value_in=(0.5, 0.9999999),
        # filter_eig_angle_in=(-180.0, -0.0001),
        filter_eig_vector_above=0.02,
        filter_eig_vector_only_in_direction=(60, 120),
        filter_eig_vector_direction_hitrate=0.1,
    )

    # Plot the selected eigenmode plane wave components (= eigen vector elements)
    if AskUserinputRecursively.yes_or_no("Plot the eigenmode components on polar?"):
        fig = MyPlotlyFigure(layout=dict(
            title_text="Absolute value of eigenmode plane wave components",
            xaxis_title="incident angle of PW [deg]",
        ))
        fig.matlab_styling()
        fig.update_layout(margin=dict(r=20, t=80, b=10))
        fig.update_polars(radialaxis=dict(range=[0, 0.1]))
        for index in indices:
            fig.add_scatterpolar(
                theta=np.rad2deg(pw_set.directions),
                r=np.abs(eig_vectors[:, index]),
                mode="markers",
                text="#%d" % index,
                name="eigenmode #%d (eigenvalue: %.3f ∠ %.3f°)" % (
                    index,
                    np.abs(eig_values[index]),
                    np.rad2deg(np.angle(eig_values[index]))
                ),
                hoverinfo="r+theta+text"
            )
        fig.show()

    # Simulate with selected eigenmode on a 2D space
    if AskUserinputRecursively.yes_or_no("Run MATLAB simulation to check the filtered eigenmodes?"):
        matlab = MatlabRunner()
        user_indices = [int(item) for item in input(
            "\nEnter eigenmode indices if you want to simulate other than the following: %s\n" % indices
        ).split()]
        if len(user_indices) != 0:
            indices = user_indices
        for index in indices:
            savemat(
                file_name="MATLAB/data/eigenvector_sim_config.mat",
                mdict=dict(
                    with_structure=True,
                    eigenmode_number=index,
                    eigen_vector=eig_vectors[:, index],
                    eval_x=np.linspace(-4.5, 4.5, 60 + 1),
                    eval_y=np.linspace(-5.0, 1.0, 40 + 1),
                )
            )
            print(
                "Running MATLAB simulation to excite with eigenmode #%d (eigenvalue: %.3f ∠ %.3f°)..." % (
                    index,
                    np.abs(eig_values[index]),
                    np.rad2deg(np.angle(eig_values[index]))
                )
            )
            matlab.run_matlab_script("mieScatt_eigenPW.m")
            matlab.export_fixer("output")
