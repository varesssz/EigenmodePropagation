import numpy as np
import pickle
from scipy.io import savemat

import AskUserinputRecursively
from MatlabRunner import MatlabRunner
from MyPlotlyFigure import MyPlotlyFigure

if __name__ == '__main__':

    # pw_set_model = "plane_wave"
    # pw_set_model = "points_along_line"
    # pw_set_model = "points_along_circle"
    pw_set_model = "points_as_phased_array"

    print(">> Initializing plane wave set from file...")
    with open("data/pw_set_from_%s.pkl" % pw_set_model, "rb") as file:
        pw_set = pickle.load(file)

    # Create indices of the excitation vector from an array of incident angles
    excitation_angles = np.deg2rad([5])[np.newaxis]
    pw_indices = np.argmin(
        np.abs(pw_set.directions - excitation_angles.T),
        axis=1,
    )
    evaluation_points_x = np.linspace(-4.5, 4.5, 240 + 1)
    evaluation_points_y = np.linspace(-5.0, 1.0, 160 + 1)

    # Save simulation configuration for MATLAB simulation
    matlab_result_fname = "output/model_%s.txt" % pw_set.model
    savemat(
        file_name="MATLAB/data/config_somePW.mat",
        mdict=dict(
            with_structure=False,
            pw_set_model=pw_set_model,
            pw_indices=pw_indices,
            eval_x=evaluation_points_x,
            eval_y=evaluation_points_y,
            result_saving_path=matlab_result_fname,
        )
    )

    if AskUserinputRecursively.yes_or_no(
            "Run MATLAB simulation to excite with partial of the PW set?\n"
            "Incident angles [deg]: %s" % [np.around(k, 2) for k in np.rad2deg(pw_set.directions[pw_indices])]
    ):
        print(">> Simulating with the model %s" % pw_set.model)
        matlab = MatlabRunner()
        matlab.run_matlab_script("mieScatt_somePW.m")
        matlab.export_fixer(matlab_result_fname)

    # Read MATLAB results
    e_z_simulated = np.genfromtxt(
        fname="MATLAB/" + matlab_result_fname,
        dtype=np.complex64,
        delimiter=",",
    )
    plotting_values = np.real(e_z_simulated)

    if AskUserinputRecursively.yes_or_no("Plot the results on a heatmap?"):
        # Create figure with heatmap plot and structure drawing
        fig = MyPlotlyFigure()
        fig.matlab_styling()
        fig.scattering_field_styling()
        fig.add_heatmap(
            x=evaluation_points_x,
            y=evaluation_points_y,
            z=plotting_values,
            colorbar_title_text="E_z(x,y) absolute value [V/m]",
            colorbar_title_side="right",
            colorscale="jet",
            zmid=0,
            # zmax=np.size(excitation_angles, 0),
            # zmin=-np.size(excitation_angles, 0),
        )
        fig.show()
