# WATER_EOS
This project combines several equations of state (EoS) of pure water, valid in localised temperature-pressure regions, to obtain one overarching EoS, spanning a very wide range in temperature and pressure.
Here is a list of the Eos combined in this work:

- [Seafreeze](https://github.com/Bjournaux/SeaFreeze/tree/master/Python), valid  in the 0-2300 MPa and 220-500 K range
- IAPWS 95, Wagner and Pruß (2002). This EoS is valid for fluid (vapor, liquid, and supercritical) phases of water, from the melting line and up to 1273 K and 1 GPa.
- Mazevet et al. (2019), valid  up to 50 000 K and hundreds of megabars
- for the region with pressure higher than the 2300 MPa we compute our own grid, based on the experimental results of Fei et al. (1993) and the theoretical effort of Mazevet  et al. (2019)

This program takes as input the temperature (in K) and pressure (in Pa) of pure water and provides several thermodynamic properties of pure water .

## Usage:
To use this EoS, import it in a python script and then initiate a DataPoint object:

```python
import EOS

my_point = EOS.DataPoint(my_temperature, my_pressure)
```
This call will generate an object containing the following list of properties:

```python
my_point.T # my_temperature, provided as input
my_point.p # my_pressure, provided as input
my_point.rho # density, mol.m^{-3}
my_point.Cp  # heat capacity at constant pressure J.mol^{-1}.K^{-1}
my_point.s  # entropy J.mol^{-1}.K^{-1}
my_point.u  # internal energy J.mol^{-1}
```
On the saturation vapor pressure line, this program returns the properties of the vapor phase. On melting lines, properties of liquid water will be provided.

The call of the EOS is further illustrated in the test_EOS.py routine.
