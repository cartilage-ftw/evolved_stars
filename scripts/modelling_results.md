# Results from modelling with spherical shell approximation

**Important summary**
* Number densities derived from the SiO $v=2$ line are at least 2 orders of magnitude lower than those from SiO $v=7$. Also, the $^{29}$ SiO $v=5, J=8-7$ transition requires $N_{mol}$ similar to SiO $v=7$. This might be a coincidence, as these are different isotopologues.
* One might need to consider how to treat maser components in this type of modelling.

Parameters that were kept constant
```python
Rstar_ep1_mas = 26.9 # in milliarcsecs, Vlemmings+17
Rstar_ep2_mas = 24.2
Tstar_ep1 = 2495
Tstar_ep2 = 2680
dist = 102 # AU
```

### Note:
I honestly wonder if there's a degeneracy between temperature, radius of the shell, and Nmol. $R_{shell}$ definitely (from experience with the model) seems to break degeneracy between $T_{exc}$, and $N_{mol}$

## SiO $v=2$
*Best fit parameters* for:
### Epoch 1
`Texc=870K, Nmol=1.35E10` for `Texc=870K, Nmol=1.35E10`
```python
deltaVelAbs = 5.8 # in km/s
deltaVelEm = 4.2 # in km/s
vstar = 39.6
absorptionShift = -9
```

#### Note
For SiO $v=2$ in epoch 1, the blueshifted (outflow) absorption component was matched. This means the $\sigma_{\rm{vel,abs}} can't be directly compared with that of epoch 2.
### Epoch 2
`Texc=700 K, Nmol=9.5E10` for an assumed `R_out=2.48*Rstar_ep2`

```python
vstar = 40.5
absorptionShift = -7.5
deltaVelAbs = 5.2 # in km/s
deltaVelEm = 5.6
```

## CO $v=1$, $J=3-2$
#### Note
For this transition's Einstein A coefficient, there were two values: `1.448e-6` (CDMS) and `2.433e-6` (SLAIM). I did the modelling with the CDMS value

### Epoch 1
Converged at `Texc=1500 K, Nmol=5E12` for `Rout = 2.066*Rstar_ep2`
with
```python
deltaVelAbs = 1.9 # in km/s
deltaVelEm = 3. # in km/s
vstar = 39.3
absorptionShift = -8.75
```

### Epoch 2
Converged at `Texc=1050 K, Nmol=15E12` for radius kept the same (`2.066*Rstar_ep2`)
```python
deltaVelAbs = 3.3 # in km/s
deltaVelEm = 3.5
vstar = 38.2
absorptionShift = -11.2
```

## SiO $v=7$
From the dirty image of epoch 1, it wasn't possible to infer a physically motivated 'radius' of the shell. I assumed `38 mas` for both.
### Epoch 1
Determined `Texc=1920 K, Nmol=1.15E12`
```python
deltaVelAbs = 1.6 # in km/s
deltaVelEm = 3.2 # in km/s
#Rstar = 2.0 #in au
vstar = 41.
absorptionShift = -7.2
```

### Epoch 2
Determined `Texc=1360 K, Nmol=3E13`
```python
deltaVelAbs = 3.1 # in km/s
deltaVelEm = 4.2 # in km/s
#Rstar = 2.0 #in au
vstar = 40.5
absorptionShift = -8.
```
## $^{29}$ SiO $v=5$
The partition function for this were taken from CDMS. Once again, I had to resort to assuming an identical radius (40 mas) for both
### Epoch 1
Determined parameters: `Texc=1350 K, Nmol=4.2E11`

**Note**: Due to the blend in the red with the H$_2$O line, it was difficult to disentangle the absorption component of the 29SiO line alone. Perhaps the velocity dispersion and frequency shift for the absorption here should be taken with a grain of salt.
```python
deltaVelAbs = 3.2 # in km/s
deltaVelEm = 4.2 # in km/s
#Rstar = 2.0 #in au
vstar = 41.7
absorptionShift = -6.
```
### Epoch 2
This line too, seemed as if there were two components. I first tried matching only one of them with `Texc=900 K, Nmol=5E12` and
```python
deltaVelAbs = 4.3 # in km/s
deltaVelEm = 3.2 # in km/s
#Rstar = 2.0 #in au
vstar = 41.
absorptionShift = -9.9
```

And then I tried to model the broader emission, arriving at `Texc=1550 K, Nmol=5E11` (notice the factor of 10 *lower* number density)
```python
deltaVelAbs = 4.1 # in km/s
deltaVelEm = 3.5 # in km/s
vstar = 38.5
absorptionShift = -12.9
```