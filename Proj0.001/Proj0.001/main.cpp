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
	int input_x;
	int input_y;
	int input_z;
	cout << "Input the number of unit cells in x,y and z direction" <<endl;
	cin >> input_x;
	cin >> input_y;
	cin >> input_z;
	Simulation* simulation2 = new Simulation(input_x,input_y,input_z,0,0,a,a,a,a,a,a,"fcc",true);
	simulation2->run_simulation();
	
	system("pause");
	return 0;
}

/*
Origo in the botom left corner, we are in the first octant! all atoms have positive coordinates. (when not moving)





*/