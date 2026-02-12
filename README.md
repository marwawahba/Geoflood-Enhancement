# Geoflood-Enhancement
The primary purpose of this supplementary material is to describe how the original GeoFlood framework was enhanced through the integration of the CPOP approach.
While the overall structure of the GeoFlood framework was preserved, several targeted adjustments were implemented to enable the application of the CPOP method.
Both the original and modified workflows are implemented in Python.
The original GeoFlood framework was developed by Zheng et. al., and the source code and accompanying guidance materials are publicly available at: https://github.com/passaH2O/GeoFlood.
Zheng et al. recommended using an Anaconda environment to run the GeoFlood modules; however, any Python-compatible environment may be used, including execution via the Windows command prompt.
Required dependencies, such as GRASS GIS and TauDEM, can be installed either prior to setting up the Python environment or configured within it.
The GitHub repository outlines the original GeoFlood workflow as a sequence of twenty steps. The modified CPOP-GeoFlood framework follows the same overall structure and numbering.
All files and steps 1 to 12 remain unchanged from the original implementation, while modifications were introduced from Step 13 onward, and the corresponding new files are provided as supplementary data.
The main objective of this document is to clearly identify and explain these code-level changes, enabling straightforward replication of the methodology presented in the main manuscript.
In addition, some modifications were introduced to support uncertainty analysis conducted prior to identifying CPOP as the preferred segmentation method.
Not all presented modifications are strictly required for implementing the CPOP-GeoFlood framework; however, they are included here for completeness and to provide guidance for users interested in extending or adapting the methodology.
