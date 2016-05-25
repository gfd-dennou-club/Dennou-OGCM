Dennou-OGCM project (Dennou Oceanic General Circulation Models project)

=======================================================================================

Target
-----------------------------------------------------------------------------------------


Models
-----------------------------------------------------------------------------------------

* Global oceanic general circulation model on aquaplanet with a spectral method (dsogcm)
  - Ocean component
    - Dynamical core: hydrostatic Boussinesq equations with a spectral Eulerian method
    - Parameterizations for sub-grid scale processes
      - Meso-scale eddy mixing
        - isopycnal diffusion scheme (Redi, 1982) 
        - GM scheme based on skew flux (Gent and McWilliams, 1990; Griffies et al, 1998)
      - Convective mixing
        - convective adjustment scheme (Marotzke, 1991)
  - Sea ice
    - Thermodynamic 3 layers sea ice model based on Winton(2000).
  - Note
    - This ocean-sea ice model can couple to an atmospheric general circulation model (DCPAM).
      For details, please see DCPCM project.
    - Coast lines can not be considered in the model
    
* Experimental models
  - Global shallow water model with a finite volume method (globalSWM-FVM)
    - The grid sysyem is triangle - hexagonal grid, and the governing equations are
      discritized with TRiSK proposed by Ringler et al. (2010).

  - Global shallow water model with a discontineous Galerkin method (globalSWM-DG)
    - The model mesh is spherical triangle mesh, and the govern equation are
      discritized with a DG method proposed by Läuter et al. (2008).


Main developers
-----------------------------------------------------------------------------------------

  - Yuta Kawai (Kobe university and RIKEN AICS)


Some informations
-----------------------------------------------------------------------------------------

- About installation
  - Please see INSTALL in the top directory.

- About license
  - Please see COPYRIGHT in the top directory.
  

Some links associated with this project
-----------------------------------------------------------------------------------------

- GFD-Dennou Club
  - https://www.gfd-dennou.org
- DCPAM: Dennou-Club Planetary Atmospheric Model
  - https://www.gfd-dennou.org/arch/dcpam/index.htm.en
- DCPCM: Dennou Club Planetary Climate simulation Model
  - https://github.com/ywkawai/DCPCM


Acknowledgments
-----------------------------------------------------------------------------------------

- Y, Kawai is receiving partial support from the Junior Research Associate (JRA) program on RIKEN. 
