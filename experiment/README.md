# Experimental Verification

## Summary
This directory contains code and data to support the experimental verification of the novel four-element LC baluns.

## Dependencies
The dependencies for the notebooks can be installed with by running `pip install -r requirements.txt` in this directory.

## Contents
- `raw_data` contains measurements of the fabricated networks as two-port Touchstone files representing each pair of ports.
- `combine_results.ipynb` is a Jupyter notebook for combining sets of two-port Touchstone files in `raw_data` to three-port Touchstone files in `processed_data`.
- `processed_data` contains measurements of the fabricated networks after being combined into three-port parameters.
- `ads_data` contains Touchstone files containing differential two-port results from Keysight Advanced Design System, either from performing electromagnetic co-simulations on the PCB designs or from simulating the Touchstone files in `processed_data`. These files have had de-embedding applied, but not renormalisation.
- `generate_figures.ipynb` performs renormalisation on the Touchstone files in `ads_data`, calculates power wave reflection coefficients and insertion loss, and displays these as figures. These figures are also saved to the `images` directory.
- `images` contains generated figures to display the simulation and experiment results.

## License
This folder and its contents are part of lc-power-match-baluns.
Copyright © 2023 Technical University of Denmark (developed by Rasmus Jepsen)

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA