import math

def get_peak_width_gradient(N, t_0, k, G):
    """Return peak width gradient case."""
    W = 4 * G * (((t_0 * (1 + k)) / math.sqrt(N)))
    return(W)


def get_peak_width_isocratic(N, t_0, k):
    """Return peak width isocratic case."""
    W = 4 * (1/math.sqrt(N)) * t_0 * (1 + k)
    return(W)


def get_G(cum_sum_w, k_phi_r, k_phi_init, t_D, t_0, t_init):
    """Return peak compression factor G."""
    G2 = (k_phi_r**2/(t_0 * (1 + k_phi_r)**2)) * ( ( ((t_D + t_init)*((1 + k_phi_init)**2))/k_phi_init**3 ) + cum_sum_w )
    G = math.sqrt(G2)
    return(G)


def poppe(k0, S, B, t0):
    """Poppe equation for peak compression factor under 1-segment linear gradient elution."""
    p = (k0*B*S*t0)/(1 + k0)
    G = math.sqrt(p**2 / 3 + p + 1)/(1+p)
    return(G)
