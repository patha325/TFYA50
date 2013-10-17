#include "simulation.h"
#include "cell_list.h"

void Simulation::run_simulation(){
	//Loop through next_time_step
	// How often shall we update cell list?

}
Simulation::Simulation (int unit_cells_x, int unit_cells_y, int unit_cells_z, int time_step,int steps,float temperature,float cutoff,float mass,float sigma,float epsilon,float lattice_constant,std::string crystal_structure,bool thermostat){
	// Create a simulation object. Initialize cell_list, list of atoms.
	create_list_of_atoms();
	create_cell_list();
	
	// Save all the input!
	
}
Simulation::~Simulation (){

}
void Simulation::create_list_of_atoms(){
	// Input
	// None
	//Requires 
	//Lattice constant & crystal_structure & unit_cells_i
	
	//Create all atoms in this class? Convert from fcc to atom positions.
	// Call atom and put in vector list_of_atoms.


	//Output
	// None

}
void Simulation::next_time_step(int current_time_step){
	/* Alter everything in the simulation to get to the next time step.
	Calculate distances
	Calculate potential
	Calculate force
	Calculate velocities
	update_atoms()
	
	{ Not every time step
	Clear cells/cell_list
	Fill cells/cell_list
	}

	Save ... to txt.

	*/
}
void Simulation::regulate_thermostat(){
	//Regulate the kinetic energy so that the temperature remains "constant"


}
void Simulation::save(){
	//Save ??? to a .txt file with some structure.

}
	

void Simulation::create_cell_list(){
	// Call a Cell_list. 
	cell_list = new Cell_list(cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant);
	cell_list->add_atoms_to_cells(list_of_atoms);
}