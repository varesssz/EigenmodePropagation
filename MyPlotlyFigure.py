from plotly.graph_objects import Figure


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
