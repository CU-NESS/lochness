# lochness
## Lunar Observatory Code in Healpy by the NESS team.

This code generates simulations of the low-frequency radio sky and celestial bodies from a given landing site on the Moon.
Users input landing site coordinate and desired time of observation in UTC and LOCHNESS will output 
the altitude/azimuth of common celestial bodies as seen from the landing site. In addition, users
can supply a healpy map to LOCHNESS and it will be rotated into the landing site reference frame
for the given UTC time. 

## Citation
If you use the algorithms or code contained in LOCHNESS, please cite (Hibbard et al?)

## GETTING STARTED:
First, make sure that all relevant dependencies are installed properly:
## Dependencies
You will need the following Python packages:
* [numpy](http://www.numpy.org/)
* [healpy](https://github.com/healpy/healpy)
* [spicepy](https://spiceypy.readthedocs.io/en/main/installation.html)

To clone a copy of the repository:
```
git clone https://github.com/CU-NESS/lochness.git
```
Then install the LOCHNESS package via:
```
cd lochness
python setup.py develop
```
Note that it is necessary to have the appropriate SPICE kernels downloaded 
(see the spicepy documentation above), including:
```
naif0012.tls,
de440.bsp,
pck00011.tpc,
moon_080317.tf,
moon_assoc_pa.tf,
moon_pa_de421_1900-2050.bpc.
```
However, for users in a hurry,
we have uploaded the most relevant ones for the LOCHNESS calculations
in the directory:
```
lochness/input/spice_kernels/
```
so they can be used immediately upon installation of LOCHNESS,
which will download the necessary files.

These are accurate and updated as of 2023.

## Contributors
Primary Authors: Joshua J. Hibbard/Valerie Wong
Secondary Authors: Neil Bassett, Keith Tauscher, Alex Eid.
