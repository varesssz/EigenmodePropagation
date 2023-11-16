import numpy as np
from scipy.io import savemat

import AskUserinputRecursively
from MatlabRunner import MatlabRunner
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

    """ Set up the simulation model """
    # Plane wave model is used if no setup is uncommented from below

    # pw_set.set_up_points_along_line_model(
    #     array_distance=-2.0 - 3.0,
    # )

    # pw_set.set_up_points_along_circle_model(
    #     array_radius=-2.0 - 3.0,
    # )

    # pw_set.set_up_phased_array_model(
    #     array_distance=-2.0 - 2.0,
    #     array_length=18,
    #     element_dist_per_lambda=0.5,
    #     taylor_windowing=True,
    # )

    # Save PW set's parameters for MATLAB usage
    pw_set.save_parameters_for_matlab(path="./MATLAB/data/")

    # Create indices of the excitation vector from an array of incident angles
    excitation_angles = np.deg2rad([120, 130])
    indices = np.argmin(
        np.abs(pw_set.directions - excitation_angles[:, np.newaxis]),
        axis=1,
    )

    # Save simulation configuration for MATLAB simulation
    matlab_result_fname = "output/model_%s.txt" % pw_set.model
    savemat(
        file_name="MATLAB/data/somePW_sim_config.mat",
        mdict=dict(
            with_structure=False,
            pw_indices=indices,
            eval_x=np.linspace(-4.5, 4.5, 240 + 1),
            eval_y=np.linspace(-5.0, 1.0, 160 + 1),
            result_saving_path=matlab_result_fname,
        )
    )

    if AskUserinputRecursively.yes_or_no(
            "Run MATLAB simulation to excite with partial of the PW set?\n%s deg" % np.rad2deg(excitation_angles)
    ):
        print("Simulating with the model %s" % pw_set.model)
        matlab = MatlabRunner()
        matlab.run_matlab_script("mieScatt_somePW.m")
        matlab.export_fixer("output")

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
        fig.add_heatmap(
            x=np.linspace(-4.5, 4.5, 240 + 1),
            y=np.linspace(-5.0, 1.0, 160 + 1),
            z=plotting_values,
            colorbar_title_text="E_z(x,y) absolute value [V/m]",
            colorbar_title_side="right",
            colorscale="jet",
            zmid=0,
            # zmax=np.size(excitation_angles, 0),
            # zmin=-np.size(excitation_angles, 0),
        )
        fig.update_layout(
            # title_text=rmsd_string,
            # margin_t=50,
            xaxis_title="x [m]",
            yaxis_title="y [m]",
            yaxis_scaleanchor="x",
            yaxis_scaleratio=1,
        )
        fig.show()
