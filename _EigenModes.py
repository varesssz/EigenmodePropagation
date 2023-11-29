import numpy as np
import pickle
from scipy.io import savemat

import AskUserinputRecursively
from MatlabRunner import MatlabRunner
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

    # pw_set_model = "plane_wave"
    # pw_set_model = "points_along_line"
    # pw_set_model = "points_along_circle"
    pw_set_model = "points_as_phased_array"

    print(">> Initializing plane wave set from file...")
    with open("data/pw_set_from_%s.pkl" % pw_set_model, "rb") as file:
        pw_set = pickle.load(file)

    cylinder_structure_name = "structureC"

    print(">> Loading Transfer-matrix from file...")
    with open("data/transfer_mat_%s_from_%s.pkl" % (cylinder_structure_name, pw_set.model), "rb") as file:
        transfer_matrix = pickle.load(file)

    print(">> Calculating eigenmodes...")
    # Eigen values and vectors of the transfer matrix
    eig_values, eig_vectors = np.linalg.eig(transfer_matrix)

    # Sort eigen values and vectors by eigen value
    sort = np.flip(np.argsort(np.abs(eig_values)))
    eig_values = eig_values[sort]
    eig_vectors = eig_vectors[:, sort]

    # Simulate with all eigenmode on a 1D space
    if AskUserinputRecursively.yes_or_no("Run MATLAB simulation to check every found eigenmode's 1D result?"):
        matlab = MatlabRunner()
        # Save simulation configuration for MATLAB simulation
        savemat(
            file_name="MATLAB/data/config_ALLeigenPW.mat",
            mdict=dict(
                pw_set_model=pw_set_model,
                eig_vectors=eig_vectors,
                eval_x=np.linspace(-4.5, 4.5, 241),
                eval_y=-0.75,
                result_saving_path="output/all_eigenmodes_by_%s" % pw_set.model,  # without file extension
            )
        )
        print(">> Running MATLAB simulation to excite with every found eigenmode...")
        print(">> Simulating with the model %s" % pw_set.model)
        matlab.run_matlab_script("mieScatt_ALLeigenPW.m")
        matlab.export_fixer("output/all_eigenmodes_by_%s.txt" % pw_set.model)

    # Select eigenmodes from the eigenvectors by filtering eigenvalues and eigenvector elements
    print(">> Filtering eigenmodes...")
    # indices = [0, 3, 6, 7, 8, 12, 18, 19, 35]
    indices = [int(item) for item in input("Enter eigenmode indices you are interested in\n(int int ...): ").split()]

    # Filer eigenvector components (alias rows) by incident angle
    direction_mask = (pw_set.directions > np.deg2rad(0)) & (np.deg2rad(180) > pw_set.directions)

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
                r=np.abs(eig_vectors[:, index]) * direction_mask,
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
        with_structure = AskUserinputRecursively.yes_or_no("Simulate with the structure in the simulation space?")
        matlab = MatlabRunner()
        user_indices = [int(item) for item in input(
            "\nEnter eigenmode indices if you want to simulate other than the following: %s\n(int int ...): " % indices
        ).split()]
        if len(user_indices) != 0:
            indices = user_indices
        for index in indices:
            savemat(
                file_name="MATLAB/data/config_eigenPW.mat",
                mdict=dict(
                    with_structure=with_structure,
                    pw_set_model=pw_set_model,
                    eigenmode_number=index,
                    eigen_vector=eig_vectors[:, index] * direction_mask,
                    eval_x=np.linspace(-4.5, 4.5, 60 + 1),
                    eval_y=np.linspace(-5.0, 1.0, 40 + 1),
                )
            )
            print(
                ">> Running MATLAB simulation to excite with eigenmode #%d (eigenvalue: %.3f ∠ %.3f°)..." % (
                    index,
                    np.abs(eig_values[index]),
                    np.rad2deg(np.angle(eig_values[index]))
                )
            )
            print(">> Simulating with the model %s" % pw_set.model)
            matlab.run_matlab_script("mieScatt_eigenPW.m")
        matlab.export_folder_fixer("output")
