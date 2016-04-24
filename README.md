# Global sea-level rise code for Ruckert et al. (in prep.)

README file last updated by Kelsey Ruckert, klr324-at-psu-dot-edu, Sun April 24 15:43:51 EST 2016

## Citation

This code is intended to accompany the results of

	Ruckert, KL, Guan, Y, Forest, FE, and Keller, K. Improving the statistical
	method can raise the upper tail of sea-level projections, (submitted to ERL).

Please cite that paper and the Rahmstorf (2007) study when using any results generated with this code.

	Rahmstorf S (2007) A Semi-empirical approach to projecting future sea-level rise
	Science 315(5810) 368–370, doi:10.1126/science.1135456.

## Overview

This code requires R with the following libraries:
- mcmc
- ncdf
- coda
- mvtnorm
- DEoptim
- compiler

This R code is intended to help users who wish to work with the sea-level rise projections or methods shown in Ruckert et al. (in prep.) in greater detail than provided in the appendix of the text. Key functionality these scripts include:

1. Global sea-level rise projections from 1880 to 2300 with associated probabilities
2. How to fit a model to observations with AR1 residuals using Markov Chain Monte Carlo and Bootstrap
3. Produces plots from the paper

The RFILES directory contains all the scripts and data necessary to run the analysis along with a README file. The Workspace subdirectory includes the prerun analysis output that was used generate the Ruckert et al. (in prep.) figures. (Note that the folder directory MUST be in the same format as when downloaded otherwise the scripts will not locate the files/scripts needed to run.)

The most important functions are **Rar** and **Robs_likelihood_AR** for MCMC calibration, **Deoptim_rahm_model** and **min_res_bootstrapped** for bootstrap calibration, and **sealevel_rahm_model**, which is the Rahmstorf (2007) sea-level model.

To fit the sea-level model to all methods, simply open R and source **Mega_Rahmstorf** (Note, this run ~5hrs). Then source **PlotRuckert_etal** to generate plots.

All of the **Toy** scripts are meant to show/teach how to fit a model using MCMC or the bootstrap method. Additionally, they were used to test bias in the surprise index. However, the program cannot be simply sourced. It requires the user to designate the length of observations and the assumptions.

----

    Copyright (C) 201X by Kelsey L. Ruckert

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
