import numpy as np
import plotly.graph_objects as go


if __name__ == '__main__':

    """ Original and subsampled E-field in total field domain """

    data = np.genfromtxt("data/sample.csv", complex, delimiter=",", skip_header=1).T
    fs = 2 * 1e9 / 3e8  # wave number sampling rate [1/m]
    window = 200  # "x" range [-100, 100] in original data [m]
    N = (int(fs * window)+1) // 2 * 2 + 1  # round up with int()+1, then make sure it's odd number
    x_less = np.linspace(-window/2, window/2, N)

    # +++ Interpolate on Descartes-coordinates +++
    # data_less_real = np.interp(x_less, np.real(data[0]), np.real(data[2]))
    # data_less_imag = np.interp(x_less, np.real(data[0]), np.imag(data[2]))
    # data_less = np.array([x_less, data[1, 0] * np.ones(N), data_less_real + data_less_imag * 1j])

    # +++ Interpolate on Polar-coordinates +++
    data_less_abs = np.interp(x_less, np.real(data[0]), np.abs(data[2]))
    data_less_ph = np.interp(x_less, np.real(data[0]), np.angle(data[2]))
    data_less = np.array([x_less, data[1, 0] * np.ones(N), data_less_abs * np.exp(1j*data_less_ph)])

    fig_Etotal = go.Figure(layout=go.Layout(title="Total E-field at sampled line"))
    fig_Etotal.update_layout(xaxis_title="x [m]", yaxis_title="E_z real value [V/m]")
    fig_Etotal.add_trace(go.Scatter(x=np.real(data[0]),
                                    y=np.real(data[2]),
                                    name="original sampling (y = 2)"))
    fig_Etotal.add_trace(go.Scatter(x=np.real(data_less[0]),
                                    y=np.real(data_less[2]),
                                    name="sub-sampled data (y = 2)"))
    # fig_Etotal.show()

    """ Plane Wave Decomposition with DFT """

    fig_fft = go.Figure(layout=go.Layout(title="Fourier-transform result"))
    fig_fft.update_layout(xaxis_title="component", yaxis_title="amplitude")
    fig_fft.add_trace(go.Scatter(y=np.abs(np.fft.fft(data[2])) / 10001,
                                 mode='markers', marker={'size': 5},
                                 name="original sampling"))

    # +++ E_z(y, k_x) with FFT{E_z(y, x)} +++
    ez_y_kx = np.fft.fft(data_less[2])
    fig_fft.add_trace(go.Scatter(y=np.abs(ez_y_kx) / N,
                                 mode='markers', marker={'size': 5},
                                 name="sub-sampled wave's decomposition"))
    fig_fft.show()

    """ Trying to reconstruct total E-field from PWD """

    k = 2 * np.pi * fs
    # +++ Creating double-sided k_x spectrum +++
    k_x = np.append(np.linspace(0, k, num=N//2+1), np.linspace(-k, 0, num=N//2+1))
    k_x = np.delete(k_x, -1)  # delete k_x=0 from the end
    k_y = np.sqrt(k**2 - k_x**2)
    # +++ E_z(y, k_x) = E_z+(k_x) * e^(+j k_y y) +++
    ez_kx_forward = ez_y_kx * np.exp(-1j * k_y * data_less[1])

    y_to_evaluate = 2
    # +++ E_z(x, y) from FFT{E_z+(k_x) * e^(+j k_y y)} +++
    ez_reconstructed = np.fft.ifft(ez_kx_forward * np.exp(1j * k_y * y_to_evaluate))
    fig_Etotal.add_trace(go.Scatter(x=np.real(data_less[0]),
                                    y=np.real(ez_reconstructed),
                                    name=f"reconstructed at (y = {y_to_evaluate})"))

    """ Validation with original simulation results from other (y = const) sampling """

    validation = np.genfromtxt(f"data/sample_validation_y0.csv", complex, delimiter=",", skip_header=1).T
    valid_N = (int(fs * 19.9)+1) // 2 * 2 + 1
    valid_x = x_less[N//2-valid_N//2:N//2+valid_N//2+1:1]  # using the same "x" values as reconstructed E_z(x, y)

    # +++ Interpolate validation data on same "x" for comparison +++
    valid_less_abs = np.interp(valid_x, np.real(validation[0]), np.abs(validation[2]))
    valid_less_ph = np.interp(valid_x, np.real(validation[0]), np.angle(validation[2]))
    valid_less = np.array([valid_x, validation[1, 0] * np.ones(valid_N), valid_less_abs * np.exp(1j*valid_less_ph)])

    fig_Etotal.update_xaxes(range=[-10, 10])  # Validation data accessible in range [-10, 10]
    fig_Etotal.add_trace(go.Scatter(x=np.real(valid_less[0]),
                                    y=np.real(valid_less[2]),
                                    name=f"validating data (y = {np.real(validation[1, 0])})"))

    fig_Etotal.show()
