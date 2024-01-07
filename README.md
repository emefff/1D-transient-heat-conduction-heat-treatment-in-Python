# 1D-transient-heat-conduction-in-Python

**This work is based on https://github.com/alexlib/trans_heat_cond/blob/master/LICENSE, thanks a lot. 
The attached code is still work in progress, it shows some weird behaviour under certain conditions (feedback on temperature by oscillating data that feeds back to data....).
To-Dos: functionalizing, cleaning up, find reason for weird oscillations etc.
Use at your own risk. Thank you!**

The 1D transient heat equation is solved for conduction and heat transfer to the surface of a steel slab. This is useful for heat treatment, especially for estimating the soaking time. It is is mostly estimated if FEA is not available (rule of thumb for example for austenitizing: 1 hour per inch of steel thickness). 
While these rules of thumb are for sure enough for a first estimate, it is not economically viable and in certain steels it can even lead to detrimental results (for example: grain coarsening in tool steels that need very high Ta).

In the shared script, a slab with a thickness of 250mm is heated in a surrounding medium of 870°C with a soaking time of 8 hours. The needed 850°C core temperature is reached already after 4 hours. Depending on how much time is needed at Ta, the soaking time may be reduced considerably.
Subsequently, the slab is cooled with different heat transfer cofficients until is reaches room temperature. Some of the input data is still not physically correct, only assumed. 

As a result we get the following temperature curves for core and surface temperature (filtered results due to oscillations in the original data):

![temperature_curves_heating_and cooling](https://github.com/emefff/Transient-heat-conduction-in-Python/assets/89903493/801444b8-4582-48cc-82e5-faf609e6b830)

From that we can easily derive the heating and cooling rates (unit: °C / minute) which we can use for comparison with TTT-diagrams etc.:

![heating_cooling_rates](https://github.com/emefff/Transient-heat-conduction-in-Python/assets/89903493/01018de6-a30a-43fc-86fb-51c9840a7e04)

emefff@gmx.at


