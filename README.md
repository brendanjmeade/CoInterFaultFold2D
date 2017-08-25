# Some notes on kinematically consistent coseismic faulting and folding in 2D

The geometric relationship between faults in folds have been extensively documented and modeled (e.g., Suppe, 1985).  The question of when do these folds grow has recieved a bit less attention (e.g., Souter and Hager, 1997; Dolan et al., 2003).  Essentially the question is whether hanging wall folds grow during earthquakes (coseismically) or in between earthquakes (interseismically).  Essentially when does deformation occur along the fold (red dashed line) in the figure below?

![Figure 1](https://github.com/brendanjmeade/CoInterFaultFold2D/blob/master/Figure1.png)

Geodetic evidence for either case appears limited but serveral models co- and inter-seismic surface deformation near thrust faults been proposed [Including but not limited to: Figure 3](https://github.com/brendanjmeade/CoInterFaultFold2D/blob/master/Figure3.png).

So if data are lacking it's the perfect time to understand the behavior of idealized fault and fold models and make predictions for the spatial patterns of surface deformation consisent with both co- and inter-seismic fold growth.  The simple models here do this for single dipping thrust fault that ramps down to a flat detachment.  A single fold bisects these two structures.  Following Souter and Hager (1997) we approximate the fold as a dislocation plane exploiting the focal mechanism ambiguity [Figure 5](https://github.com/brendanjmeade/CoInterFaultFold2D/blob/master/Figure5.png).  This means that deformation associated with material passing through the kink axis is treated elastically.

The individual constributions to a surface deformation are shown in ![Figure 6](https://github.com/brendanjmeade/CoInterFaultFold2D/blob/master/Figure6.png) with blue lines representing horizontal motion and green lines representing vertical motion.  The upper right panel shows the long term motions over many earthquake cycles with steps at the surface trace of both the fault and fold.  The region between the fault and fold in uniformly uplifted.  The upper right panel shows surface motions in response to unit slip on the fault ramp.  The lower left panel shows the surface motions in response to unit slip on a small portion of the horizontal detachment.  The lower right panel shows the surface motions in response to unit slip across the fold axis.  Over geologic time the last three panels must sum to match the first.

The predicted surface displacements from the different classes of models are shown in: ![Figure 9](https://github.com/brendanjmeade/CoInterFaultFold2D/blob/master/Figure9.png)

The 6 panels show: a) interseimsic horizontal deformation (note that the label is incorrect in this panel!), b) interseismic vertical deformation, c) coseismic horizontal deformation, d) coseismic vertical deformation, e) cumulative horizontal deformation, f) cumulative vertical deformation. Predictions from the four different kinematic models (definined in figure 3) colored as 1-green, 2-blue, 3-cyan, 4-red.

The central results may be summarized as:

1 - The clearest distinctions between the different kinematic models can be seen in the botton two panels (e, f) which show the total cumulative (or geologic) deformation.  The kinematically consistent fault and fold model (red line) produces distinct displacement jumps at *both* the locations where both the fault and fold intersect the surface rather than just at the fault trace (other 3 models).  In particular for the flat detachment case considered here long-term vertical uplift is localized entirely *between* the surface intersection between the fault and fold and is consistent with long-term antiformal growth. 

2 - Difference in coseismic displacements are localized around the intersection of the surface with the fold axis.

All figures (except figures 1 and 5) are produced by the matlab scripts included here.
