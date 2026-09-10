# MW/M31 Halo Major Merger Remnant Structure 

Density profile of the merged Dark Matter remnant following the predicted Milky Way and M31 major merger.

ASTR 400B research project, University of Arizona, January to May 2025. The report is written in LaTeX using MNRAS formatting.


---
## Contents

| File | File overview | 
|---|---|
| FinalCode_revised.ipynb | Revised analysis. Corrected center of mass, binning, and R200. |
| FinalCodE.ipynb | Original submission, unedited | 
| TORRES_MASTR400B_Final.pdf | Written report, formatted in MNRAS style | 
| proj.tex | LaTeX source | 
| NFW1.py | Navarro-Frenk-White density profile | 
| Hernquist1.py | Hernquist density profile | 
| CenterOfMass.py, <br> CenterOfMass2.py | Center of mass classes. Version 2 uses the shrinking sphere method. | 
| ReadFile.py, Read1.py | Snapshot file readers | 
| Lab7-Copy1.ipynb | 2D contour plot of the merged remnant, showing how the dark matter particles settle after the major merger | 
| MW_000.txt, M31_000.txt, <br> M31_779.txt | Simulation snapshot data |

Note: MW_779.txt exceeds GitHub's file size limit so it is not included.

---
## Objective
Analyze the final remnant structural evolution of the Dark Matter halo during a major merger between 
two galaxies, comparing the merged halo's density profile to Navarro-Frenk-White and Hernquist Profiles
to investigate how accurately these analytic models represent Dark Matter behavior.

---
## Methods

N-body simulation data (Milky Way–M31, snapshot 779, ~11.1 Gyr), Hernquist and NFW analytical models, 
Center of Mass computed with the shrinking sphere method for all three halos, radial density binning, 
R200 computation from the total density contained within each radius. Python.

---
## Results

- The NFW profile closely tracks the simulated density distribution, holding even at large radii
   with a median deviation of 13%, compared to a median ratio of 3.07 for Hernquist.
- Computed a well defined remnant halo boundary of R200 = 242.5 kpc.
- Inner regions remained dense and stable while outer regions experienced the greatest redistribution.
- A 2D contour map of the remnant shows the halo densest at the core with a symmetrical particle distribution,
 and density decreasing at larger radii.

---
## Outcome

The NFW profile is a superior fit to the simulated density distribution after the collision, supporting 
the CDM prediction that Dark Matter halos follow a universal density profile despite experiencing a major
merger. The inner region remains stable while the outer region undergoes tremendous redistribution, consistent 
with gravity keeping the core condensed while the outer halo is pushed outwards. It also demonstrates how structural 
changes during mergers provide the foundation for decoding Dark Matter's true nature.

One limitation is the choice of radial binning, which can distort how the redistribution of particles is evaluated. 
Another is the selection of scale radius, which can shift the density profile and restrict how closely the models can be tested. 
The analytic models also expect specific conditions to match the halos but ignore physical forces such as dynamical friction. 
Binning and scale radius were held constant across both models so the scale radius remained the dependent variable. 

---
## Note on revision

Revised September 2026. A single weighted average had pulled the merged remnant's center away from the
density peak, flattening the inner profile. It now uses the shrinking sphere method, as the individual 
halos already did. All profiles use a single radial bin array so the comparison between the merged halo 
and the individual galaxies is consistent. R200 is computed from the total density contained within each 
radius instead of the density in each shell, which corrected it from 150.8 to 242.5 kpc.

The report was written before these corrections and uses the original R200 value.
Its conclusion identifies binning, scale radius, and a shifted center of mass as the main limitations, 
which is where these corrections came from.
