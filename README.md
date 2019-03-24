DCPOM project (Dennou-Club Planetary Oceanic Model project)

=======================================================================================

Target
-----------------------------------------------------------------------------------------
DCPOM project designs and develops ocean and sea ice models in order to explore climates on exoplanets.

Models
-----------------------------------------------------------------------------------------

* Global oceanic general circulation model for ocean planet climates (dogcm)
  - Ocean component
    - Dynamical core: hydrostatic Boussinesq equations
    - Parameterizations for sub-grid scale processes
      - Meso-scale eddy mixing
        - isopycnal diffusion scheme (Redi, 1982) 
        - GM scheme based on skew flux (Gent and McWilliams, 1990; Griffies et al, 1998)
      - Convective mixing
        - convective adjustment scheme (Marotzke, 1991)
    - Numerical method
      - horizontal discretization: a spectral Eulerian method based spherical harmonics 
      - vertical discretization: a finite volume method
      - time discretization
        - a LF-AM3 scheme (Shchepetkin and McWilliams, 2005) for advective process
	- a leapfrog scheme for fast wave process and vertical diffusive process,
	  a backward Euler scheme for lateral diffusive process
  - Sea ice component
    - Thermodynamic 3 layers sea ice model based on a formulation by Winton(2000).
  - Note
    - This ocean-sea ice model can couple to an atmospheric general circulation model (DCPAM).
      For details, please see DCPCM project.
    - Coast lines can not be considered in the model
    
* Experimental models
  - Global shallow water model with a finite volume method (globalSWM-FVM)
    - The grid sysyem is triangle - hexagonal grid, and the governing equations are discretized with TRiSK proposed by Ringler et al. (2010).

  - Global shallow water model with a discontinuous Galerkin method (globalSWM-DG)
    - The model mesh is spherical triangle mesh, and the governing equation are discretized with a DG method proposed by Läuter et al. (2008).


Some information
-----------------------------------------------------------------------------------------

- About installation
  - Please see INSTALL in the top directory.

- About license
  - Please see COPYRIGHT in the top directory.
  
- About what is contained in each directory
  - Please see DIR_LIST in the top directory. 


Main developers
-----------------------------------------------------------------------------------------

  - Yuta Kawai (Kobe university and RIKEN R-CCS) and the collaborators in GFD-Dennou Club


Some links associated with this project
-----------------------------------------------------------------------------------------

- GFD-Dennou Club
  - https://www.gfd-dennou.org
- DCPAM: Dennou-Club Planetary Atmospheric Model
  - https://www.gfd-dennou.org/arch/dcpam/index.htm.en
- DCPCM: Dennou-Club Planetary Climate Model
  - https://github.com/gfd-dennou-club/Dennou-CCM
- DCPOM-DG: Experimental OGCM with a discontinuous galerkin method 
  - https://github.com/ywkawai/DCPOM-DG

Acknowledgments
-----------------------------------------------------------------------------------------

- Yuta Kawai received partial support from the Junior Research Associate (JRA) program on RIKEN (FY2015-FY2017). 
