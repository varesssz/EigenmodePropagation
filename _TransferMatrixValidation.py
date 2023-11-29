import numpy as np
import pickle
from scipy.io import savemat

import AskUserinputRecursively
from MatlabRunner import MatlabRunner
from MyPlotlyFigure import MyPlotlyFigure
from PlaneWaveSet import PlaneWaveSet

if __name__ == '__main__':

    print(">> Initialization...")
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

    print(">> Loading Transfer-matrix from file...")
    with open("data/transfer_mat_%s_from_%s.pkl" % (cylinder_structure_name, pw_set.model), "rb") as file:
        transfer_matrix = pickle.load(file)

    # Create excitation vector
    excitation_angles = np.deg2rad([10, 15, 20, 80, 90, 91, 140, 180])
    pw_indices = np.argmin(
        np.abs(pw_set.directions - excitation_angles[:, np.newaxis]),
        axis=1,
    )
    excitation_vector = np.zeros(pw_set.k_x.size, np.complex64)
    excitation_vector[pw_indices] = 1 + 0j

    # Calculate Ez(k_x) with transfer-matrix and excitation vector
    print(">> Calculating resulted plane waves from transfer-matrix...")
    ez_kx = transfer_matrix @ excitation_vector

    # Reconstruct total E-field with this result
    y_reconstruct = 0.0
    print(">> Reconstructing total E-field from plane wave components at y=%.1fm..." % y_reconstruct)
    ez_reconstructed_1d = np.zeros(pw_set.x_sampling.size, dtype=np.complex64)
    for k_x_i, k_y_i, ez_kx_forward_i in zip(pw_set.k_x, pw_set.k_y, ez_kx):
        ez_reconstructed_1d += (
            ez_kx_forward_i * np.exp(1j * k_y_i * y_reconstruct) * np.exp(1j * k_x_i * pw_set.x_sampling)
        )

    # Save simulation configuration for MATLAB simulation
    savemat(
        file_name="MATLAB/data/somePW_sim_config.mat",
        mdict=dict(
            with_structure=True,
            pw_indices=pw_indices,
            eval_x=pw_set.x_sampling,
            eval_y=y_reconstruct,
            result_saving_path="data/validation_e_field_y0.txt",
        )
    )

    if AskUserinputRecursively.yes_or_no(
            "Run MATLAB simulation to excite with partial of the PW set?\n"
            "%s deg" % [np.around(k, 2) for k in np.rad2deg(pw_set.directions[pw_indices])]
    ):
        print(">> Simulating with the model %s" % pw_set.model)
        matlab = MatlabRunner()
        matlab.run_matlab_script("mieScatt_somePW.m")
        matlab.export_fixer("data/validation_e_field_y0.txt")

    print(">> Reading MATLAB results...")
    e_z_simulated = np.genfromtxt(
        fname="MATLAB/data/validation_e_field_y0.txt",
        dtype=np.complex64,
        delimiter=",",
    )

    # Root-Mean-Square Deviation
    print(">> Calculating RMS Deviation of simulated and reconstructed E-field real values...")
    rmsd = (np.real(e_z_simulated) - np.real(ez_reconstructed_1d))**2
    rmsd = np.sqrt(rmsd.sum() / rmsd.size)
    nrmsd = rmsd / (np.real(e_z_simulated).max() - np.real(e_z_simulated).min())
    rmsd_string = "RMSD: %s,     NRMSD: %.2f%%  (%s)" % (
        np.format_float_scientific(rmsd, 2), nrmsd * 100, np.format_float_scientific(nrmsd, 2)
    )
    print(rmsd_string)

    # Adding reconstruction and validation to total E-field plot for comparison
    print(">> Creating Plotly figure with the reconstructed and the validation total E-field plots...")
    fig = MyPlotlyFigure()
    fig.update_layout(
        xaxis_title="x [m]",
        yaxis_title="E_z real value [V/m]",
    )
    fig.add_scatter(
        x=pw_set.x_sampling,
        y=np.real(ez_reconstructed_1d),
        name=f"reconstructed",
        mode="markers",
    )
    fig.add_scatter(
        x=pw_set.x_sampling,
        y=np.real(e_z_simulated),
        name=f"simulated",
    )
    fig.matlab_styling()
    fig.update_layout(title_text=rmsd_string, margin_t=50)
    fig.update_xaxes(range=[-10, 10])
    fig.show()
