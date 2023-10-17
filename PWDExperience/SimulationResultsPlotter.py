import numpy as np
from MyPlotlyFigure import MyPlotlyFigure


if __name__ == '__main__':
    data2 = np.genfromtxt(
        "data/single_cylinder_E_field_10400x10400mm@20x20mm_f1000MHz.txt",
        skip_header=2,
    ).T

    # data = np.array([
    #     data2[0, 187560:188081:]/1000,
    #     data2[1, 187560:188081:]/1000,
    #     data2[7, 187560:188081:] + 1j * data2[8, 187560:188081:],
    # ])

    data2_heatmap = data2[7] + 1j * data2[8]
    data2_heatmap = data2_heatmap.reshape((521, 521))
    data2_heatmap_x = data2[0].reshape((521, 521))
    data2_heatmap_y = data2[1].reshape((521, 521))
    data2_fig = MyPlotlyFigure()
    data2_fig.matlab_styling()
    data2_fig.update_layout(
        xaxis_title="x [m]",
        yaxis_title="y [m]",
        # plot_bgcolor='white',
        # font_size=18,
        # margin=dict(r=20, t=20, b=10),
    )
    data2_fig.add_heatmap(
        z=np.abs(data2_heatmap)[110:410:, 60:461:],
        x=np.real(data2_heatmap_x)[0, 60:461:] / 1000,
        y=np.real(data2_heatmap_y)[110:410:, 0] / 1000,
        colorscale="rainbow",
        colorbar_title="E_z absolute value [V/m]",
        colorbar_titleside="right",
        zmin=0,
        # zmid=0,
        zmax=1.75,
    )
    data2_fig.add_scatter(
        x=np.array([0]),
        y=np.array([0]),
        mode="markers",
        marker=dict(
            size=60, color="rgba(255, 255, 255, 0)",
            line=dict(
                color="rgba(255, 255, 255, 1)",
                width=2
            ),
        ),
    )
    data2_fig.show()
