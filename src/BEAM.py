import numpy as np
import matplotlib.pyplot as plt
import os 
wave_lgh = 633e-9              # He-Ne (m)
w0 = 1e-3                      # waist (m)
L = 10e-3                      # Size env (m)
N = 1024                       # Number of point


def gaussian_beam(n , size , w0):
    x = np.linspace(-size / 2 , size/2 , n , endpoint=False)
    X , Y = np.meshgrid (x,x)
    return np.exp(-(X**2 + Y**2)/w0**2), x


def propagate(u0 , size , wave_lgh, i ):
    n = u0.shape[0]
    fx = np.fft.fftfreq(n , d = size /n)
    FX , FY = np.meshgrid(fx , fx)
    arg = 1 - (wave_lgh * FX )** 2 - (wave_lgh * FY)**2
    kz = (2 * np.pi/wave_lgh) * (np.sqrt(np.maximum(arg , 0)))
    dil_phase = np.exp(1j * kz * i) * ( arg > 0)         # Phase delay of the wave after passing through z
    return np.fft.ifft2(np.fft.fft2(u0) * dil_phase)


def beam_width(u , x):
    Intensity = (np.abs(u) **2).sum (axis = 0)
    sigma_dispersion = (Intensity * x**2).sum() / Intensity.sum()
    return 2 * np.sqrt(sigma_dispersion)


def main():
    u0 , x = gaussian_beam (N , L , w0)
    rayleigh = np. pi * w0**2 /wave_lgh
    z_point = np.linspace(0 , 2 * rayleigh , 20)

    w_num = [beam_width(propagate(u0 , L , wave_lgh , i), x) for i in z_point]
    w_z = w0 * np.sqrt(1+(z_point / rayleigh)**2)

    # FIGURE PART
    os.makedirs("figures", exist_ok=True)
    plt.figure(figsize=(6,4))
    plt.plot(z_point , np.array(w_z) * 1e3 , label = "Analytical")
    plt.plot(z_point , np.array(w_num) * 1e3 , "o" , markersize = 5 , label = "Angular spectrum")
    plt.xlabel("Z (m)")
    plt.ylabel("Beam Radius w (mm)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("figures/beam_width.png" , dpi = 200)
    fig , axes = plt.subplots(1 , 3 , figsize = (10 , 3.2))

    for ax, z in zip(axes, [0, rayleigh, 2 * rayleigh]):
        intensity = np.abs(propagate(u0 , L , wave_lgh , z ))**2
        ax.imshow(intensity, extent=[x[0] * 1e3, x[-1] * 1e3] * 2, cmap="inferno")
        ax.set_title(f"z = {z :.1f} m")
        ax.set_xlabel("x (mm)")
    plt.tight_layout()
    plt.savefig("figures/beam_profiles.png", dpi=200)
    print("Done")

if __name__ == "__main__":
    main()


