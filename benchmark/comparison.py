import os
import subprocess

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

elements = np.logspace(2, 7, 6)

time_nim = []
time_fenics = []

fenics_env = os.environ.copy()
fenics_env["OMP_NUM_THREADS"] = "1"

os.makedirs("timings", exist_ok=True)
run = True
for n in [int(np.sqrt(n / 2)) for n in elements]:
    print(n)

    if run:
        subprocess.run(["/usr/bin/time", "-o", f"timings/nim{n}.txt", "./poisson", f"{n}"])
        subprocess.run(
            [
                "/usr/bin/time",
                "-o",
                f"timings/fenics{n}.txt",
                "python3",
                "poisson.py",
                f"{n}",
            ],
            env=fenics_env,
        )

    with open(f"timings/nim{n}.txt") as f:
        line = f.readline()
        user_time = line.split("user")[0]
        time_nim.append(float(user_time))

    with open(f"timings/fenics{n}.txt") as f:
        line = f.readline()
        user_time = line.split("user")[0]
        time_fenics.append(float(user_time))

time_nim = np.array(time_nim)
time_fenics = np.array(time_fenics)
plt.loglog(elements, time_nim, marker=".", label="Nimfem")
plt.loglog(elements, time_fenics, marker=".", label="FEniCS")
plt.legend()
plt.xlabel("Number of Elements")
plt.ylabel("Runtime [s]")

plt.savefig("poisson_comparison.png")
