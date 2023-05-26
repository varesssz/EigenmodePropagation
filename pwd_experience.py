import warnings

import numpy as np
import plotly.graph_objects as go


if __name__ == '__main__':

    """ Plotting Figures """
    debug_plotter = go.Figure()  # for debug reasons, its data gets flushed before plotting
    fig_fft = go.Figure()  # for FFT results plotting, ex. compare before and after sampling data
    fig_Etotal = go.Figure()  # for the total E-filed plotting, ex. compare original and reconstructed field

    """ ----------------------------------------------------- """
    """ Original and subsampled E-field in total field domain """
    """ ----------------------------------------------------- """

    """ Reading samples from file """
    data = np.genfromtxt(
        "data/sample_y2.csv",
        complex,
        delimiter=",",
        skip_header=1,
    ).T

    """ Parameters of PWD """
    wavelength = 299792458 / 1e9  # [m]
    k_wave = 2 * np.pi / wavelength
    sampling_rate = 14 / wavelength  # [1/m]
    window = 120  # "x" range [-100, 100] in original data [m]
    if np.abs(data[0, 0] - data[0, 1]) > 1 / sampling_rate:
        warnings.warn("Oversampling the dataset!", RuntimeWarning)

    """ Numbers of the samples rounded up with int()+1, then making sure it's odd number """
    N = (int(window * sampling_rate) + 1) // 2 * 2 + 1

    """ sampling is made to continuous wave likeness """
    x_sampling = np.linspace(-window/2, window/2, N, endpoint=False)
    x_step = np.mean(np.delete(np.delete(np.append(x_sampling, 0) - np.append(0, x_sampling), -1), 0))
    if x_step > wavelength / 2:
        warnings.warn("Undersampling the expected wavelength!", RuntimeWarning)

    """ Interpolate on Descartes-coordinates """
    # data_sampled_real = np.interp(x_sampling, np.real(data[0]), np.real(data[2]))
    # data_sampled_imag = np.interp(x_sampling, np.real(data[0]), np.imag(data[2]))
    # data_sampled = np.array(
    #     [x_sampling,
    #      data[1, 0] * np.ones(N),
    #      data_sampled_real + data_sampled_imag * 1j]
    # )

    """ Interpolate on Polar-coordinates """
    data_sampled_abs = np.interp(x_sampling, np.real(data[0]), np.abs(data[2]))
    data_sampled_ph = np.interp(x_sampling, np.real(data[0]), np.angle(data[2]))
    data_sampled = np.array(
        [x_sampling,
         data[1, 0] * np.ones(N),
         data_sampled_abs * np.exp(1j * data_sampled_ph)]
    )

    """ Plotting total E-field """
    fig_Etotal.update_layout(
        # title="Total E-field at sampled line",
        xaxis_title="x [m]",
        yaxis_title="E_z real value [V/m]",
    )
    fig_Etotal.add_trace(
        go.Scatter(
            x=np.real(data[0]),
            y=np.real(data[2]),
            name="original sampling (y = 2)",
        )
    )
    fig_Etotal.add_trace(
        go.Scatter(
            x=np.real(data_sampled[0]),
            y=np.real(data_sampled[2]),
            name="sub-sampled data (y = 2)",
        )
    )
    fig_Etotal.show()

    """ --------------------------------- """
    """ Plane Wave Decomposition with DFT """
    """ --------------------------------- """

    """ E_z(k_x, y) with FFT{E_z(x, y)}"""
    ez_kx_y = np.fft.fft(data_sampled[2]) / N

    """ Plotting FFT results """
    fig_fft.update_layout(
        title="Fourier-transform result",
        xaxis_title="component",
        yaxis_title="amplitude"
    )
    fig_fft.add_trace(
        go.Scatter(
            y=np.real(np.fft.fft(data[2])) / 10001,
            mode='markers', marker={'size': 5},
            name="original sampling",
        )
    )
    fig_fft.add_trace(
        go.Scatter(
            y=np.real(ez_kx_y),
            mode='markers', marker={'size': 5},
            name="sub-sampled wave's decomposition",
        )
    )
    fig_fft.show()

    """ k_x_max = 2pi/∆x * (N-1)/N as maximum of FFT frequency"""
    k_x_step = (2 * np.pi / x_step) * (1 / N)
    """ Creating positive k_x[i] array from 0 to k_x_max/2 and negative k_x[i] array from -k_x_max/2 to 0 """
    k_x = np.concatenate((
        np.linspace(
            0,
            k_x_step * (N-1) / 2,
            int(N / 2 + 1),  # rounding up half of odd N
        ),
        np.linspace(
            -k_x_step * (N-1) / 2,
            0,
            int(N / 2),  # rounding down half of odd N
            endpoint=False,
        )
    ))
    """ delete redundant k_x=0 from the array's end """
    # k_x = np.delete(k_x, -1)  # used without "endpoint=False" flag at the last np.linspace()

    """ Throwing away the decaying plane wave components """
    """ Also sorting in ascending order """
    ez_kx_y = np.append(
        ez_kx_y[int(N - k_wave / k_x_step + 1):],
        ez_kx_y[:int(k_wave / k_x_step + 1)]
    )
    k_x = np.append(
        k_x[int(N - k_wave / k_x_step + 1):],
        k_x[:int(k_wave / k_x_step + 1)]
    )
    debug_plotter.data = []
    debug_plotter.add_trace(go.Scatter(y=np.real(ez_kx_y), mode="markers"))
    debug_plotter.add_trace(go.Scatter(y=np.imag(ez_kx_y), mode="markers"))
    debug_plotter.show()

    """ Creating k_y from k_x """
    k_y = np.sqrt(k_wave**2 - k_x**2)

    """ Attempting to use decaying PW components too """
    """ Probably just needed to get it in order (neg,0,pos instead of 0,pos,neg) """
    # k_y_pos_prop = np.sqrt(k_wave**2 - k_x[(k_x < k_wave) & (k_x >= 0)]**2)
    # k_y_decay = -1j * np.sqrt(-k_wave**2 + k_x[np.abs(k_x) > k_wave]**2)
    # k_y_neg_prop = np.sqrt(k_wave**2 - k_x[(k_x > -k_wave) & (k_x < 0)]**2)
    # k_y = np.concatenate((k_y_pos_prop, k_y_decay, k_y_neg_prop))

    debug_plotter.data = []
    debug_plotter.add_trace(go.Scatter(y=np.real(k_x), mode="markers"))
    debug_plotter.add_trace(go.Scatter(y=np.real(k_y), mode="markers"))
    debug_plotter.add_trace(go.Scatter(y=np.imag(k_y), mode="markers"))
    debug_plotter.show()

    """ Compensating that the sampling is not starting from x=0 coordinate"""
    ez_kx_y_xcomp = ez_kx_y * np.exp(1j * k_x * data_sampled[0, 0])
    debug_plotter.data = []
    debug_plotter.add_trace(go.Scatter(y=np.real(ez_kx_y_xcomp), mode="markers"))
    debug_plotter.add_trace(go.Scatter(y=np.imag(ez_kx_y_xcomp), mode="markers"))
    debug_plotter.show()

    """ Compensate that the sampling is not in the y=0 line """
    """ with E_z(k_x,y) = E_z_forward(k_x) * exp(+j k_x x) as E_z_backward can be negligible in E_z(k_x,y) """
    ez_kx_forward = ez_kx_y_xcomp * np.exp(-1j * k_y * data_sampled[1, 0])
    debug_plotter.data = []
    debug_plotter.add_trace(go.Scatter(y=np.real(ez_kx_forward), mode="markers"))
    debug_plotter.add_trace(go.Scatter(y=np.imag(ez_kx_forward), mode="markers"))
    debug_plotter.show()

    """ ----------------------------------------------------- """
    """ Reconstruction in 1 dim of the total E-field from PWD """
    """ ----------------------------------------------------- """

    x_range = [-15, 15]
    y_reconstruct = 2.0

    x_reconstruct = x_sampling[
        int(N / 2 + N * x_range[0] / window):
        -int(N / 2 - N * x_range[1] / window) + 1:
    ]

    ez_reconstructed = np.zeros(np.shape(x_reconstruct), dtype=np.complex_)
    for k_x_i, k_y_i, ez_kx_forward_i in zip(k_x, k_y, ez_kx_forward):
        ez_reconstructed = ez_reconstructed + \
                           ez_kx_forward_i * np.exp(1j * k_y_i * y_reconstruct) * np.exp(1j * k_x_i * x_reconstruct)

    """ Adding to total E-field plot for comparison """
    fig_Etotal.add_trace(
        go.Scatter(
            x=np.real(x_reconstruct),
            y=np.real(ez_reconstructed),
            name=f"reconstructed (y = {y_reconstruct})",
        )
    )
    fig_Etotal.update_xaxes(range=x_range)
    fig_Etotal.show()

    """ ----------------------------------------------------- """
    """ Reconstruction in 2 dim of the total E-field from PWD """
    """ ----------------------------------------------------- """

    x_range = [-6, 6]
    y_range = [-3, 3]

    y_reconstruct = np.linspace(
        y_range[0],
        y_range[1],
        num=int((y_range[1]-y_range[0]) / x_step),
        endpoint=False
    )
    y_reconstruct = y_reconstruct[:, np.newaxis]
    x_reconstruct = x_sampling[
        int(N / 2 + N * x_range[0] / window):
        -int(N / 2 - N * x_range[1] / window) + 1:
    ]

    ez_reconstructed = np.zeros(np.shape(x_reconstruct), dtype=np.complex_)
    for k_x_i, k_y_i, ez_kx_forward_i in zip(k_x, k_y, ez_kx_forward):
        ez_reconstructed = ez_reconstructed + \
                           ez_kx_forward_i * np.exp(1j * k_y_i * y_reconstruct) * np.exp(1j * k_x_i * x_reconstruct)

    heatmap = go.Figure()
    heatmap.add_trace(
        go.Heatmap(
            x=np.real(x_reconstruct),
            y=np.real(y_reconstruct)[:, 0],
            z=np.real(ez_reconstructed),
            colorscale="rainbow",
            zmid=0,
        )
    )
    heatmap.show()

    """ --------------------------------------------------------------------------- """
    """ Validation with original simulation results from other (y = const) sampling """
    """ --------------------------------------------------------------------------- """

    validation = np.genfromtxt(
        f"data/sample_y2.csv",
        complex,
        delimiter=",",
        skip_header=1
    ).T

    valid_window = 40
    valid_N = (int(valid_window * sampling_rate) + 1) // 2 * 2 + 1

    """ Using the same "x" values as reconstructed E_z(x, y) """
    valid_x = x_sampling[N // 2 - valid_N // 2:N // 2 + valid_N // 2 + 1:1]

    """ Interpolate validation data on same "x" for comparison"""
    valid_sampled_abs = np.interp(valid_x, np.real(validation[0]), np.abs(validation[2]))
    valid_sampled_ph = np.interp(valid_x, np.real(validation[0]), np.angle(validation[2]))
    valid_samping = np.array(
        [valid_x,
         validation[1, 0] * np.ones(valid_N),
         valid_sampled_abs * np.exp(1j * valid_sampled_ph)]
    )

    """ Adding validation data to total E-field plot """
    fig_Etotal.add_trace(
        go.Scatter(
            x=np.real(valid_samping[0]),
            y=np.real(valid_samping[2]),
            name=f"validating data (y = {np.real(validation[1, 0])})",
        )
    )

    """ Set plot appearance for report/thesis figures """
    fig_Etotal.update_layout(
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

    fig_Etotal.show()
