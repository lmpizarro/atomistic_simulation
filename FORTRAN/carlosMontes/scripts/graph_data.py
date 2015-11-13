import numpy as np
import matplotlib.pyplot as plt

N = 500
a = np.loadtxt("presion.030.txt")

print np.mean(a[N:]), np.var(a[N:]), np.std(a[N:])

plt.plot(a[N:])

plt.show()


hist, bins = np.histogram(a[N:], bins=50)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width=width)
plt.show()
