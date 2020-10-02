import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# consts
N = 500

# before the quench
mu0 = 0
t0 = 1
d0 = -1

# after the quench
mu1 = -6
t1 = t0
d1 = d0

# nonlinearity strength
lmd = 0

# small time number
R = 160

# large time number
Q = 1000

# small time step
dt = 0.00025

# large time step

ds = R * dt

tol = 1e-16

cutOff = 1.2

timeAxisParts=10

dirVal = "/home/disk2/Documents/cppCode/kerrKitaev/benchmark/"
datName = "mu0" + str(mu0) + "t0" + str(t0) + "d0" + str(d0) + "mu1" + str(mu1) + "t1" + str(
    t1) + "d1" + str(d1) + "lmd" + str(lmd)
inGFileName = dirVal + "G" + datName + ".csv"

inMat = pd.read_csv(inGFileName, header=None)

inMat = inMat.iloc[range(0, int(N / 2)), :]
rowN = len(inMat)
colN = len(inMat.columns)

outMat = np.zeros((rowN, colN))
for i in range(0, rowN):
    for j in range(0, colN):
        outMat[i, j] = inMat.iloc[i, j] % (np.pi)

fig, ax = plt.subplots(1, 1)

img = ax.imshow(outMat, extent=[0, 1, 0, 1])

xInd = np.arange(0, Q + Q / timeAxisParts, Q / timeAxisParts) * ds

x_label_list = xInd
xTickList = [ind / (Q * ds) for ind in xInd]
ax.set_xticks(xTickList)

ax.set_xticklabels(x_label_list)

yIndRange = np.arange(0, timeAxisParts+1)
yMax = max(yIndRange)
yTickList = [n / yMax for n in yIndRange]
ax.set_yticks(yTickList)
y_label_list = [str(n / yMax) + "$\pi$" for n in yIndRange]
ax.set_yticklabels(y_label_list)

fig.colorbar(img)
plt.title("Geometric phase, $\mu_{0}=$" + str(mu0) + ", $t_{0}=$" + str(t0) + ", $\Delta_{0}=$" + str(d0) + ", $\mu_{1}=$" + str(
    mu1) + ", $t_{1}=$" + str(t1) + ", $\Delta_{1}=$" + str(d1) + ", $\lambda=$" + str(lmd))
plt.xlabel("time")
plt.ylabel("k")

outPicName = dirVal + "phase" + datName + ".png"
plt.savefig(outPicName)
plt.close()
