import numpy as np
import math
import pandas as pd
from IPython.display import display

# Parameters
e_values = [0.0167, 0.2]
Delta_f = math.pi/5
f_targets = np.array([k*Delta_f for k in range(0, 11)])  # 0 to 2π
t0, tf = 0.0, 1.0
N = 5000
h = (tf - t0)/N

def C_of_e(e):
    return 2*math.pi*(1-e*2)*(-1.5)

def f_dot(t, f, e):
    return C_of_e(e) * (1 + e * math.cos(f))**2

def rk4_integrate(e):
    t = t0
    f = 0.0
    ts = [t]
    fs = [f]
    for _ in range(N):
        k1 = f_dot(t, f, e)
        k2 = f_dot(t+h/2, f + h*k1/2, e)
        k3 = f_dot(t+h/2, f + h*k2/2, e)
        k4 = f_dot(t+h, f + h*k3, e)
        f = f + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        t = t + h
        ts.append(t)
        fs.append(f)
    return np.array(ts), np.array(fs)

# function to interpolate t at target f
def interpolate_t(ts, fs, f_vals):
    t_at_f = []
    for ft in f_vals:
        idx = np.searchsorted(fs, ft)
        if idx == 0:
            t_at_f.append(ts[0])
        elif idx >= len(fs):
            t_at_f.append(ts[-1])
        else:
            t1, t2 = ts[idx-1], ts[idx]
            f1, f2 = fs[idx-1], fs[idx]
            # linear interpolation
            t_interp = t1 + (ft - f1)*(t2 - t1)/(f2 - f1)
            t_at_f.append(t_interp)
    return t_at_f

tables = {}

for e in e_values:
    ts, fs = rk4_integrate(e)
    t_interp = interpolate_t(ts, fs, f_targets)
    df = pd.DataFrame({
        "f (rad)": f_targets,
        "t correspondente": t_interp
    })
    tables[e] = df
    display(f"Tabela e={e}", df)

tables