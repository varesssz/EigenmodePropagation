import numpy as np
from MyPlotlyFigure import MyPlotlyFigure
import AskUserinputRecursively


if __name__ == '__main__':

    # pw_set_model = "plane_wave"
    # pw_set_model = "points_along_line"
    # pw_set_model = "points_along_circle"
    pw_set_model = "points_as_phased_array"

    # Which structure to draw on the heatmap
    cylinder_structure_name = "structureC"
    cylinder_struct = np.genfromtxt(
        fname="data/cylinders_%s.txt" % cylinder_structure_name,
        dtype=np.float64,
        delimiter=",",
        skip_header=2,
    )

    # Check the following eigenmodes (plotting can be skipped on each iteration)
    indices = [int(item) for item in input("Enter eigenmode indices you are interested in\n(int int ...): ").split()]
    for index in indices:
        print("\n\n>> Eigenmode #%d" % index)
        folder_path = "MATLAB/data/structC_w120_x5605_k801/points_along_circle_model/"
        # folder_path = "MATLAB/output/"

        # Import PEC simulation results
        e_field_pec = np.genfromtxt(
            fname=folder_path + "eigenmode_by_%s_e_field_%d_pec.txt" % (pw_set_model, index),
            dtype=np.complex64,
            delimiter=",",
        )
        plotting_values = np.abs(e_field_pec)

        if AskUserinputRecursively.yes_or_no("Plot the eigenmode on heatmap?"):
            # Create figure with heatmap plot and structure drawing
            fig = MyPlotlyFigure()
            fig.matlab_styling()
            fig.scattering_field_styling(cylinder_struct=cylinder_struct)
            fig.add_heatmap(
                x=np.linspace(-4.5, 4.5, 60 + 1),
                y=np.linspace(-5.0, 1.0, 40 + 1),
                z=plotting_values,
                colorbar_title_text="E_z(x,y) absolute value [V/m]",
                colorbar_title_side="right",
            )
            fig.show()

        if AskUserinputRecursively.yes_or_no("Calculate the difference from vacuum simulation?"):
            # Import Vacuum simulation results
            e_field_vacuum = np.genfromtxt(
                fname=folder_path + "eigenmode_by_%s_e_field_%d_vacuum.txt" % (pw_set_model, index),
                dtype=np.complex64,
                delimiter=",",
            )
            plotting_values = np.abs(e_field_pec) - np.abs(e_field_vacuum)

            # Root-Mean-Square Deviation
            rmsd = (np.abs(e_field_pec) - np.abs(e_field_vacuum)) ** 2
            rmsd = np.sqrt(rmsd.sum() / rmsd.size)
            nrmsd = rmsd / (np.abs(e_field_vacuum).max() - np.abs(e_field_vacuum).min())
            rmsd_string = "RMSD: %s,     NRMSD: %.2f%%  (%s)" % (
                np.format_float_scientific(rmsd, 2), nrmsd * 100, np.format_float_scientific(nrmsd, 2)
            )
            print(rmsd_string)

            if AskUserinputRecursively.yes_or_no("Plot the difference on heatmap?"):
                # Create figure with heatmap plot and structure drawing
                fig = MyPlotlyFigure()
                fig.matlab_styling()
                fig.scattering_field_styling(cylinder_struct=cylinder_struct)
                fig.update_layout(
                    title_text=rmsd_string,
                    margin_t=50,
                )
                fig.add_heatmap(
                    x=np.linspace(-4.5, 4.5, 60 + 1),
                    y=np.linspace(-5.0, 1.0, 40 + 1),
                    z=plotting_values,
                    colorbar_title_text="E_z(x,y) absolute value difference [V/m]",
                    colorbar_title_side="right",
                    colorscale="jet",
                    zmax=np.abs(e_field_pec).max(),
                    zmin=-np.abs(e_field_pec).max(),
                )
                fig.show()
