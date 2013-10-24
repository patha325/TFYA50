#include "simulation.h"
#include <iostream>
#include "vec.h"
using namespace std;


int main(){
    
	//initializegraphics
	float a = 0.2; // Watch out for double warning!
    
/*
Parameters:
	int unit_cells_x				: Number of unit cells in x-direction
	int unit_cells_y				: Number of unit cells in y-direction
	int unit_cells_z				: Number of unit cells in x-direction
	int time_step					: Length of one time step
	int steps						: number of time steps for the simulation
    float temperature				: Simulation temperature
    float cutoff					: Cutoff distance for potential calculations
    float mass						: Mass of an atom
    float sigma						: Material constant sigma
    float epsilon					: Material constan epsilon
    float lattice_constant			: Lattice constant
    std::string crystal_structure	: Name of chrystal structure (fcc, hcp, bcc...)
    bool thermostat					: If a thermostat is employed
 */
	Simulation* simulation2 = new Simulation(8,4,4,0,0,a,4.4,a,a,a,4.3,"Hej!",true);
	simulation2->run_simulation();
	
	//system("pause");
	return 0;
}

/*
Origo in the botom left corner, we are in the first octant! all atoms have positive coordinates. (when not moving)





*/