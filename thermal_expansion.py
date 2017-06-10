import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

# Helped Article : http://stackoverflow.com/questions/28974818/heat-diffusion-on-a-2d-plate-python

def set_heatsource(M, Heat):
    # function to set M values corresponding to non-zero Gr Value
    M[0:L, 49:52] = np.where(Heat > 0, Heat, M[0:L, 49:52])


def black_body_radiation(T, area, specific_heat, dt):
    '''
    This function will return temperature change induced by black body radiation.
    :param T: Temperature of matrix material 
    :param area: Contact area between adjacent pixel that assumed cubic shape. so pow(dx,2) could be area.
    :param specific_heat: Specific heat of matrix material
    :param dt: Time interval between each step
    :return: 
    '''
    sigma = 5.67 * 0.00000001
    Wloss = (pow(T, 4) - pow(298, 4)) * area * sigma
    Tloss = (Wloss * dt) / specific_heat
    return Tloss


def temp_loss_by_air(Tmat, Tair, dx, tc, area, specific_heat, dt):
    """
    as a result of this function, temperature loss by air with thermal conduction will returned.
    :param Tmat: Temperature of Matrix
    :param Tair: Temparature of Air
    :param dx: Size of each pixel
    :param tc: Thermal Conductivity
    :param area: Contact area between adjacent pixel that assumed cubic shape. so pow(dx,2) could be area.
    :param specific_heat: Specific heat of matrix material
    :param dt: Time interval between each step
    :return:
    """
    # calculate energy loss by air. value of '0.4' in this function is ambigous value.
    energyflow = area * ((Tmat - Tair) / ((dx / tc) + (dx / 0.4))) * dt
    Tloss = energyflow / specific_heat
    return Tloss


# simulation configurations
dt = 4  # step interval of simulation, sec
dx = 0.01  # size of each side of voxel, meter
time = 20000  # Total duration of simulation running time, sec
Temp = 298  # Initial temperature of 2-dimensional plate, Kelvin

# material properties
sh = 0.4479876  # specific heat (J/g*K)
tc = 150  # thermal conductivity (W/mK)
density = 7860000  # (g/m^3)
te = 12.5 * 0.000001  # (m/mK)

# shape properties of matrix
area = dx * dx  # (m^2)
volume = dx * dx * dx  # (m^3)
mass = volume * density  # g
L = 5  # Length of the Plate
B = 101  # Width of the Plate

# shape of heating device. initial temperature set 1K for easy manipulation of temperature by multiplying numbers.
HSource = np.ones([L, 3])

# print simulation information when simulation starts
print("simulation grid size : " + str(dx * L) + "m*" + str(dx * B) + "m")
print("total simulation time : " + str(time) + "sec")
print("step interval : " + str(dt) + "sec")

plt.style.use('ggplot')

# matrix generation with given temperature '273K'
M = np.ones([L, B]) * Temp
set_heatsource(M, HSource)

# Build MM, a list of matrices, each element corresponding to M at a given step
T = np.arange(0, time, dt)  # evenly spaced timeframe from 0 to 10 with 0.1 interval
Matrix = []
LM = []
DIL = []
TIME = []
pointT = []

for i in range(len(T)):
    for j in range(0, B - 1):
        for i in range(0, L - 1):
            ##set boundary condition
            if 0 == i and 44 < j < 57:
                M[0, j] = M[1, j]
            elif L - 1 == i and 44 < j < 57:
                M[L - 1, j] = M[L - 2, j]
            else:
                M[0, j] = M[1, j] - black_body_radiation(M[L - 1, j], area, sh, dt) - temp_loss_by_air(M[L - 1, j], 298,
                                                                                                       dx, tc, area, sh,
                                                                                                       dt)
                M[L - 1, j] = M[L - 2, j] - black_body_radiation(M[L - 2, j], area, sh, dt) - temp_loss_by_air(
                    M[L - 2, j], 298, dx, tc, area, sh, dt)
                M[i, 0] = M[i, 1] - black_body_radiation(M[i, 1], area, sh, dt) - temp_loss_by_air(M[i, 1], 298, dx, tc,
                                                                                                   area, sh, dt)
                M[i, B - 1] = M[i, B - 2] - black_body_radiation(M[i, B - 2], area, sh, dt) - temp_loss_by_air(
                    M[i, B - 2], 298, dx, tc, area, sh, dt)

            ##determine the each voxel boundary temperature difference
            tT = M[i + 1, j] - M[i, j]  # top
            lT = M[i, j - 1] - M[i, j]  # left
            rT = M[i, j + 1] - M[i, j]  # right
            bT = M[i - 1, j] - M[i, j]  # bottom

            ##calculate estimated energy flow (J/s = W)
            topFlow = (tT * tc * (area / dx))
            leftFlow = (lT * tc * (area / dx))
            botFlow = (bT * tc * (area / dx))
            rightFlow = (rT * tc * (area / dx))

            ##calculate temperature change
            # interaction with left interface
            dlT = (leftFlow / (sh * mass))  # temperature difference between pixel
            M[i, j - 1] = M[i, j - 1] - dlT

            # interaction with top interface
            dtT = (topFlow / (sh * mass))  # temperature difference between pixel
            M[i + 1, j] = M[i + 1, j] - dtT

            # interaction with bottom interface
            dbT = (botFlow / (sh * mass))  # temperature difference between pixel
            M[i - 1, j] = M[i - 1, j] - dbT

            # interaction with right interface
            drT = (rightFlow / (sh * mass))  # temperature difference between pixel
            M[i, j + 1] = M[i, j + 1] - drT

            M[i, j] = M[i, j] + drT + dlT + dtT + dbT
    shape = np.shape(Matrix)

    # setting temperature profile.
    if shape[0] * dt < 3600:
        HT = Temp + (((shape[0] * dt) / (60)) * 10)
        HSource = np.ones([L, 1])
        HSource = HSource * HT  # temperature arise with 10K/min heating rate
        set_heatsource(M, HSource)
        print('Heating, T=' + str(HT))


    elif 3600 <= shape[0] * dt <= 7200:
        HSource = np.ones([L, 1])
        HT = Temp + (((3600) / (60)) * 10)
        HSource = HSource * HT  # holding temperature on specific value.
        set_heatsource(M, HSource)
        print('Holding, T=' + str(HT))

    else:
        print('Cooling')
        # no energy input, whole specimen cooled by thermal radiation and conduction.
    Matrix.append(M.copy())

    # get thermal dilatation from temperature profile
    LP = M[4]  # linear temperature profile(axial center for maximum dilatation)
    PTP = M[0, 50]  # point temperature profile(Surface)
    LK = LP - (np.ones(B) * Temp)
    TD = sum(LK * te * dx)
    DIL.append(TD.copy())  # Thermal Dilatation Profile
    LM.append(LP.copy())  # LM Temperature Profile
    TIME.append(shape[0] * dt)  # Elapsed Time from Simulation Start
    pointT.append(PTP)

    print(str(shape[0]) + "/" + str(time / dt) + " completed")

fig, ax1 = plt.subplots()
ax1.plot(TIME, DIL, 'b-')
ax1.set_xlabel('time(s)')
ax1.set_ylabel('thermic dilatation(m)', color='b')
ax2 = ax1.twinx()
ax2.plot(TIME, pointT, 'r.')
ax2.set_ylabel('temperature(K)', color='r')
ax2.set_ylim([200, 1000])
fig = plt.figure()
pcm = plt.pcolormesh(Matrix[200])
plt.colorbar()


## Function called for updating the graphic by timely manner
def step(i):
    if i >= len(Matrix): return
    pcm.set_array(Matrix[i].ravel())
    plt.draw()


anim = FuncAnimation(fig, step, interval=1)  # .save('simul.mp4', writer="ffmpeg")
plt.show()
