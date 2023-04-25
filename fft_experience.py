import numpy as np
import plotly.graph_objects as go

if __name__ == '__main__':

    # Playing with DFT and sampling in time domain

    N = 48
    t_sample = np.linspace(0, 1, N + 1)
    t_analog = np.linspace(0, 1, 100*N+1)
    w = 2*np.pi
    f = 3
    y_sample = np.cos(f * w * t_sample) + np.cos(2 * f * w * t_sample) + np.sin(2/3 * f * w * t_sample)
    y_analog = np.cos(f * w * t_analog) + np.cos(2 * f * w * t_analog) + np.sin(2/3 * f * w * t_analog)

    fig_time = go.Figure(layout=go.Layout(title="Time Domain"))
    fig_time.update_layout(xaxis_title="time [s]", yaxis_title="Amplitude")
    fig_time.add_trace(go.Scatter(x=t_analog, y=y_analog,
                                  name="Analog"))
    fig_time.add_trace(go.Scatter(x=t_sample, y=y_sample,
                                  mode='markers', marker={'size': 10},
                                  name="Samples"))

    fft = np.fft.fft(y_sample, N)

    fig_fft = go.Figure(layout=go.Layout(title="Frequency Domain"))
    fig_fft.update_layout(xaxis_title="frequency [1/s]", yaxis_title="Amplitude")
    fig_fft.add_trace(go.Scatter(y=np.real(fft)/N,
                                 mode='markers', marker={'size': 10},
                                 name="Fourier REAL"))
    fig_fft.add_trace(go.Scatter(y=np.imag(fft)/N,
                                 mode='markers', marker={'size': 10},
                                 name="Fourier IMAG"))
    fig_fft.show()

    y_re = np.fft.ifft(fft)

    fig_time.add_trace(go.Scatter(x=t_sample, y=np.real(y_re),
                                  name="Reconstructed"))
    fig_time.show()
