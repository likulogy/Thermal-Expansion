# Thermal-Expansion
This is my first step for computer assisted engineering :D 

I wrote this code for personal study on area of Computer Simulation :)

Purpose of this simulation is get thermal expansions in single materials under condition of conduction and black body radiation.

All of assumes, and functions are quite ambigious. so DO NOT TRUST result of this simulation.

# Considered Parameters
## Thermal Conduction
```python
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
```
Thermal conduction calculated by energy flow that induced by temperature difference between adjacent pixels.

First, calculate temperature difference between pixel, like `tT = M[i + 1, j] - M[i, j]  # top`.

Second, calculate energy flow induced by temperature difference based on thermal conduction equation. like `topFlow = (tT * tc * (area / dx))`.

If speciman have higher thermal conductivity 'tc', energy flow will increased. also, if contact area decreased energy flow will retarded also.

After all these procedure, calculate amount of temperature change by division with specific heat of matrix material. see below.

```python
##calculate temperature change
# interaction with left interface
dlT = (leftFlow / (sh * mass))  # temperature difference between pixel
M[i, j - 1] = M[i, j - 1] - dlT
```
## Thermal Radiation
Thermal radiation considered by really simple manner, the Stefan-Boltzmann law. I assume thermal radiation have an proportionality with fourth power of temperature.

```python
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
```
Using defined function above, I able to calculate temperature loss by black body radiation that have proportionality with Stefan-Boltzmann law. 

NOTE : Influence of area change by thermal dilatation does not considered in this simulation. but you can easily adjust simulation.

## Thermal Dilatation
Thermal Dilatation aquired from temperature line-profile.

First, I aquire line temperature profile through whole matrix(`LP = M[4]`). Then, simply calculate a deviation from 'Initial Temperature'. The deviation value is stored in 'LK' array. 

Then, simply multiplicate with Coefficient of Thermal Expansion.


```python
LK = LP - (np.ones(B) * Temp)
TD = sum(LK * te * dx)
```

finally, I could get an value of Thermal Dilatation, 'TD'.

# Result
![result demonstration](/result_sim.png?raw=true "Optional Title")

# Acknowledgement
This code got lots inspiration and help from http://stackoverflow.com/questions/28974818/heat-diffusion-on-a-2d-plate-python
