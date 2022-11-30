#space 
./launchLUMI.sh convergence.sh fluid1 0.4 1.0 0.1 5 4
./launchLUMI.sh convergence.sh fluid2 0.4 1.0 0.1 6 4
./launchLUMI.sh convergence.sh fluid3 0.4 1.0 0.1 7 4
./launchLUMI.sh convergence.sh fluid4 0.4 1.0 0.1 7 6
./launchLUMI.sh convergence.sh fluid5 0.4 1.0 0.1 8 4
./launchLUMI.sh convergence.sh fluid6 0.4 1.0 0.1 9 4

#Rtol
./launchLUMI.sh convergence.sh rtol1 0.4 0.1 0.01 7 6
./launchLUMI.sh convergence.sh rtol2 0.4 1.0 0.10 7 6
./launchLUMI.sh convergence.sh rtol3 0.4 1.5 0.15 7 6
./launchLUMI.sh convergence.sh rtol4 0.4 3.0 0.30 7 6
./launchLUMI.sh convergence.sh rtol5 0.4 4.0 0.40 7 6
./launchLUMI.sh convergence.sh rtol6 0.4 10. 1.00 7 6

#CFL
./launchLUMI.sh convergence.sh cfl1 0.1 1.0 0.1 7 6
./launchLUMI.sh convergence.sh cfl2 0.2 1.0 0.1 7 6
./launchLUMI.sh convergence.sh cfl3 0.3 1.0 0.1 7 6
./launchLUMI.sh convergence.sh cfl4 0.4 1.0 0.1 7 6
./launchLUMI.sh convergence.sh cfl5 0.5 1.0 0.1 7 6
./launchLUMI.sh convergence.sh cfl6 0.8 1.0 0.1 7 6
