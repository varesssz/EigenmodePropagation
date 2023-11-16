import numpy as np
from MyPlotlyFigure import MyPlotlyFigure
import AskUserinputRecursively


if __name__ == '__main__':

    # Which structure to draw on the heatmap
    cylinder_structure_name = "structureC"
    cylinder_struct = np.genfromtxt(
        fname="data/cylinders_%s.txt" % cylinder_structure_name,
        dtype=np.float64,
        delimiter=",",
        skip_header=2,
    )

    # Check the following eigenmodes (plotting can be skipped on each iteration)
    for index in [757, 761, 763, 768, 772, 776]:  # [757, 761, 763, 768, 772, 776]
        print("\n\n>> Eigenmode #%d" % index)
        folder_path = "MATLAB/data/structC_w120_x5605_k801/0deg_to_180deg_excitation/"
        # folder_path = "MATLAB/output/"

        # Import PEC simulation results
        e_field_pec = np.genfromtxt(
            fname=folder_path + "eigenmode_e_field_%d_pec.txt" % index,
            dtype=np.complex64,
            delimiter=",",
        )
        plotting_values = np.abs(e_field_pec)

        if AskUserinputRecursively.yes_or_no("Calculate the difference from vacuum simulation?"):
            # Import Vacuum simulation results
            e_field_vacuum = np.genfromtxt(
                fname=folder_path + "eigenmode_e_field_%d_vacuum.txt" % index,
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

        if AskUserinputRecursively.yes_or_no("Plot the eigenmode on heatmap?"):
            # Create figure with heatmap plot and structure drawing
            fig = MyPlotlyFigure()
            fig.matlab_styling()
            fig.add_heatmap(
                x=np.linspace(-4.5, 4.5, 240 + 1),
                y=np.linspace(-5.0, 1.0, 160 + 1),
                z=plotting_values,
                colorbar_title_text="E_z(x,y) absolute value [V/m]",
                colorbar_title_side="right",
                # colorscale="jet",
                # zmid=0,
                zmax=14,
                zmin=0,
            )
            fig.add_scatter(
                x=cylinder_struct[:, 0],
                y=cylinder_struct[:, 1],
                mode="markers",
                marker=dict(
                    size=45, color="rgba(255, 255, 255, 0)",
                    line=dict(
                        color="rgba(255, 255, 255, 1)",
                        width=2
                    ),
                ),
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
