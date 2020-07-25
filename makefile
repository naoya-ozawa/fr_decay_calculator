decaychain:	decaychain.cpp
	`root-config --cxx --cflags` -o decaychain decaychain.cpp `root-config --glibs`

flux_calculator:	flux_calculator.cpp
	`root-config --cxx --cflags` -o flux_calculator flux_calculator.cpp `root-config --glibs`

espectrum:	espectrum.cpp
	`root-config --cxx --cflags` -o espectrum espectrum.cpp `root-config --glibs`
