from plotly.graph_objects import Figure
from numpy import ndarray


class MyPlotlyFigure(Figure):
    def matlab_styling(self):
        self.update_layout(
            plot_bgcolor='white',
            font_size=18,
            margin=dict(r=20, t=20, b=10),
            xaxis=dict(
                showline=True,
                linecolor="grey",
                showgrid=True,
                gridcolor="lightgrey",
                zerolinewidth=2,
                zerolinecolor="lightgrey"
            ),
            yaxis=dict(
                showline=True,
                linecolor="grey",
                showgrid=True,
                gridcolor="lightgrey",
                zerolinewidth=2,
                zerolinecolor="lightgrey"
            ),
        )
        self.update_polars(
            bgcolor="white",
            sector=(0, 180),
            angularaxis=dict(
                showline=True,
                linecolor="grey",
                showgrid=True,
                gridcolor="lightgrey",
            ),
            radialaxis=dict(
                showline=True,
                linecolor="grey",
                showgrid=True,
                gridcolor="lightgrey",
            ),
        )

    def scattering_field_styling(self, cylinder_struct: ndarray = None):
        self.update_layout(
            xaxis_title="x [m]",
            yaxis_title="y [m]",
            yaxis_scaleanchor="x",
            yaxis_scaleratio=1,
        )
        if cylinder_struct is not None:
            self.add_scatter(
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
