import numpy as np
import plotly.graph_objects as go


if __name__ == '__main__':

    """ Original and subsampled E-field in total field domain """

    data = np.genfromtxt("data/sample.csv", complex, delimiter=",", skip_header=1).T
    wavelength = 299792458 / 1e9
    k_wave = 2 * np.pi / wavelength
    resolution = wavelength / 20
    window = 60  # "x" range [-100, 100] in original data [m]
    N = (int(window / resolution)+1) // 2 * 2 + 1  # round up with int()+1, then make sure it's odd number
    x_sampling = np.linspace(-window/2, window/2, N+1)
    x_sampling = np.delete(x_sampling, -1)  # Im{gs}
    x_step = x_sampling[1] - x_sampling[0]

    # +++ Interpolate on Descartes-coordinates +++
    # data_less_real = np.interp(x_less, np.real(data[0]), np.real(data[2]))
    # data_less_imag = np.interp(x_less, np.real(data[0]), np.imag(data[2]))
    # data_less = np.array([x_less, data[1, 0] * np.ones(N), data_less_real + data_less_imag * 1j])

    # +++ Interpolate on Polar-coordinates +++
    data_sampled_abs = np.interp(x_sampling, np.real(data[0]), np.abs(data[2]))
    data_sampled_ph = np.interp(x_sampling, np.real(data[0]), np.angle(data[2]))
    data_sampled = np.array([x_sampling, data[1, 0] * np.ones(N), data_sampled_abs * np.exp(1j * data_sampled_ph)])

    fig_Etotal = go.Figure(layout=go.Layout(title="Total E-field at sampled line"))
    fig_Etotal.update_layout(xaxis_title="x [m]", yaxis_title="E_z real value [V/m]")
    # fig_Etotal.add_trace(go.Scatter(x=np.real(data[0]),
    #                                 y=np.real(data[2]),
    #                                 name="original sampling (y = 2)"))
    # fig_Etotal.add_trace(go.Scatter(x=np.real(data_sampled[0]),
    #                                 y=np.real(data_sampled[2]),
    #                                 name="sub-sampled data (y = 2)"))
    # fig_Etotal.show()

    """ Plane Wave Decomposition with DFT """

    k_x = np.append(np.linspace(0, (N-1) / N * np.pi / x_step, int((N+1) / 2)),
                    np.linspace(-(N-1) / N * np.pi / x_step, -2 * np.pi / x_step / N, int((N-1) / 2)))
    k_x = np.append(k_x[round(N-k_wave/np.pi*x_step*(N-1)/2):], k_x[:round(k_wave/np.pi*x_step*(N-1)/2+1)])  # Re{ks}
    k_y = np.sqrt(k_wave**2 - k_x**2)  # Im{ks}
    # +++ E_z(y, k_x) with FFT{E_z(y, x)} +++
    ez_y_kx = np.fft.fft(data_sampled[2]) / N

    fig_fft = go.Figure(layout=go.Layout(title="Fourier-transform result"))
    fig_fft.update_layout(xaxis_title="component", yaxis_title="amplitude")
    fig_fft.add_trace(go.Scatter(y=np.real(np.fft.fft(data[2])) / 10001,
                                 mode='markers', marker={'size': 5},
                                 name="original sampling"))
    fig_fft.add_trace(go.Scatter(y=np.real(ez_y_kx),
                                 mode='markers', marker={'size': 5},
                                 name="sub-sampled wave's decomposition"))
    # fig_fft.show()

    ez_kx = np.append(ez_y_kx[round(N - k_wave / np.pi * x_step * (N - 1) / 2):],
                      ez_y_kx[: round(k_wave / np.pi * x_step * (N - 1) / 2 + 1)])

    ez_kx_forward = ez_kx * np.exp(-1j * k_y * data_sampled[1, 0]) * np.exp(1j * k_x * data_sampled[0, 0])

    """ Trying to reconstruct total E-field from PWD """

    # k = 2 * np.pi * fs/2
    # # +++ Creating double-sided k_x spectrum +++
    # k_x = np.append(np.linspace(0, k, num=N//2+1), np.linspace(k, 0, num=N//2+1))
    # k_x = np.delete(k_x, -1)  # delete k_x=0 from the end
    # k_y = np.sqrt(k**2 - k_x**2)
    # +++ E_z(y, k_x) = E_z+(k_x) * e^(+j k_y y) +++

    y_to_evaluate = 0
    # +++ E_z(x, y) from FFT{E_z+(k_x) * e^(+j k_y y)} +++
    # ez_reconstructed = np.fft.ifft(ez_kx_forward * np.exp(1j * k_y * y_to_evaluate))

    # ez_reconstructed = np.array([])
    # for i in np.arange(np.size(x_sampling)):
    #     summa = 0
    #     for n in np.arange(np.size(k_x)):
    #         summa += ez_kx_forward[n] * np.exp(1j*k_y[n]*y_to_evaluate) * np.exp(1j*k_x[n]*x_sampling[i])
    #     ez_reconstructed = np.append(ez_reconstructed, summa)

    ez_reconstructed = np.array(np.shape(x_sampling), dtype=np.complex_)
    for k_x_i, k_y_i, ez_kx_forward_i in zip(k_x, k_y, ez_kx_forward):
        ez_reconstructed = ez_reconstructed + ez_kx_forward_i * np.exp(1j*k_y_i*y_to_evaluate) * np.exp(1j*k_x_i*x_sampling)
    ez_reconstructed -= N

    fig_Etotal.add_trace(go.Scatter(x=np.real(data_sampled[0]),
                                    y=np.real(ez_reconstructed),
                                    name=f"reconstructed at (y = {y_to_evaluate})"))

    """ Validation with original simulation results from other (y = const) sampling """

    validation = np.genfromtxt(f"data/sample_validation_y0.csv", complex, delimiter=",", skip_header=1).T
    valid_N = (int(19.9 / resolution)+1) // 2 * 2 + 1
    valid_x = x_sampling[N // 2 - valid_N // 2:N // 2 + valid_N // 2 + 1:1]  # using the same "x" values as reconstructed E_z(x, y)

    # +++ Interpolate validation data on same "x" for comparison +++
    valid_less_abs = np.interp(valid_x, np.real(validation[0]), np.abs(validation[2]))
    valid_less_ph = np.interp(valid_x, np.real(validation[0]), np.angle(validation[2]))
    valid_less = np.array([valid_x, validation[1, 0] * np.ones(valid_N), valid_less_abs * np.exp(1j*valid_less_ph)])

    fig_Etotal.update_xaxes(range=[-10, 10])  # Validation data accessible in range [-10, 10]
    fig_Etotal.add_trace(go.Scatter(x=np.real(valid_less[0]),
                                    y=np.real(valid_less[2]),
                                    name=f"validating data (y = {np.real(validation[1, 0])})"))

    fig_Etotal.show()
