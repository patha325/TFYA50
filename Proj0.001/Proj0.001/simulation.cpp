#include "simulation.h"

using namespace std;

/*-----------------------
CONSTRUCTOR
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
-
Creates a simulation object. 
Calls constructors for all atoms and the cell list. 
------------------------*/
Simulation::Simulation (int unit_cells_x, int unit_cells_y, int unit_cells_z, int time_step,int steps,float temperature,float cutoff,float mass,float sigma,float epsilon,float lattice_constant,std::string crystal_structure,bool thermostat){
	create_list_of_atoms();
	create_cell_list();
	Vec test1 (1,2,3);
	Atom p (test1);
	Vec test2 (0,0,0);
	Atom m (test2);
		
	// Todo: Save all the input!	
}

/*-------------------------
DESTRUCTOR
Destroys all atoms and the 
cells and cell list.
-------------------------*/
Simulation::~Simulation (){

}

/*----------------------------
FUNCTION: run_simulation
Paramteters: None
Returns: Nothing
-
Here is the main loop for the 
simulation. Handles time steps, 
and everything that happes during
the simulation.
----------------------------*/
void Simulation::run_simulation(){

	// Todo: How often shall we update cell list?
}

/*-----------------------------
FUNCTION: create_list_of_atoms
Parameters: None
Returns: Nothing
-
Create all atoms and add them to the
vector list_of_atoms.
-----------------------------*/
void Simulation::create_list_of_atoms(){
	/*
	// Calculate number of atoms
		Vec extra (0,0,0);
	if(crystal_structure == "fcc"){
	for(int i=1;i<=unit_cells_x;i++){
		Vec origin (0,0,0);
		list_of_atoms.push_back(new Atom(origin+extra,0.5));
		extra= Vec (lattice_constant,0,0);

	}
	*/
	}
	//Todo: Lattice constant & crystal_structure & unit_cells_i
	//Todo: Create all atoms in this class? Convert from fcc to atom positions.
}

/*------------------------------
FUNCTION next_time_step
Paramteters: int
Returns: Nothing
-
Alter everything in the simulation to get to the next time step.
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
------------------------------*/
void Simulation::next_time_step(int current_time_step){}

/*------------------------------
FUNCTION regulate_thermostat
Paramteters: None
Returns: Nothing
- 
Regulates systems to maintain the 
temperature of the heat sink.
------------------------------*/
void Simulation::regulate_thermostat(){}

/*----------------------------
FUNCTION save
Parameters: None
Return: Nothing
-
Save ??? to a .txt file with some structure. 
TODO: How save, format?
----------------------------*/
void Simulation::save(){}
	
/*------------------------------
FUNCTION create_cell_list
Parameters: None
Returns: Nothing
-
Creates the cell list by calling its 
constructor. Adds all atoms to the 
cell list.
------------------------------*/
void Simulation::create_cell_list(){

	cell_list = new Cell_list(cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant);
	cell_list->add_atoms_to_cells(list_of_atoms);
}