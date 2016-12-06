# SESL: Semi-Empirical Sea Level 

README file last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Dec 05 10:26:00 EST 2016

## Citation and Acknowledgements

Version 1.0 of this code was released on 5 December 2016 to accompany

	Kopp, R. E., A. C. Kemp, K. Bittermann, B. P. Horton,
	J. P. Donnelly, W. R. Gehrels, C. C. Hay, J. X. Mitrovica,
	E. D. Morrow, and S. Rahmstorf (2016). Temperature-driven
	global sea-level variability in the Common Era. Proceedings
	of the National Academy of Sciences. doi: 10.1073/pnas.1517056113.
	
Please cite this source when using this code.

This code with developed by Klaus Bittermann, with assistance from Eric Morrow.

## Folders

- Data: holds all input data
- Simu: holds the output of 'Calc_SL_from_Param.m'
- Out: holds all other output

## Main Code

- Run_SESL.m: Runs everything and gives the opportunity to define all settings

- SESL.m: Main script with MH-algorithm which outputs the posterior parameter set
- DefineSettings_SESL.m: Defines all settings for SESL.m
- LoadData_SESL.m: Loads the input data for SESL.m

- Calc_SL_from_Param.m: Calculates sl, T0, c & temperature time series from the posterior parameter set (SESL.m output) and saves them. The number of files to split the output to can be set in order to save memory.

- Calc_SESL_Prc.m: Calculates percentiles of sea level and other parameters necessary for projections. To do this it loads the files stored by Calc_SL_from_Param.m above.

- Calc_SESLProjection.m: Calculates projections with RCP temperature inputs from the output of Calc_SESL_Prc.m

- Calc_SESLConterfact.m: Calculates sea level for counterfactual temperature inputs from the output of Calc_SESL_Prc.m

Other external functions called by SESL.m & Calc_SL_from_Param.m:

- calc_sl
- calc_temp
- calc_T0
- calc_prob
- find_optim_Ho
- resize_T


----

    Copyright (C) 2016 by Klaus Bittermann

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
