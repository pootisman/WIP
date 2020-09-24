# WIParse - generate pretty plots from Wireless Insite SBR results!

## Description

A small toolkit for analyzing/plotting the simulation results from Wireless Insite software. Apart from that, a set of utilities for shifting storage onto online DB is added for convenience. Package is still under development (contributions are welcome!) and is planned to be a subsystem of another tool for performing radio studies and simulations (not the actual simulation engine!).

## Capabilities

* Multithreaded SQLite database uploader for MySQL databases (with authentication)
* CIR constructor and collator for later comparison
* Logarithmic model coefficient estimation (WIP)
* Received signal pattern polar plot (2D only, 3D is planned)
* Channel image plot
* RX/TX group filtering as well as capability of setting indexes of receivers
* All data for SISO is stored (MIMO planned)

## TODO

* 3D mesh storage (Grab model references and floor plans from project?)
* Temporal/Spatial characteristics analysis