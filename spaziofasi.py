import matplotlib.pyplot as plt
import numpy as np

x, v = np.loadtxt("dati.dat", unpack=True, usecols=(1, 2))
plt.plot(x, v, color="g", linestyle="-")
plt.title("Spazio delle fasi")
plt.xlabel("Posizione");
plt.ylabel("Velocità");
plt.savefig("ultimospaziofasi.png")
plt.show()
