#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>

#include "atom.h"
#include "cell_list.h"
#include "cell_list.h"
#include "vec.h"

class Simulation{
private:
	//Variables
	std::vector<Atom*> list_of_atoms;
	int number_of_atoms;
	float time_step; //Given in femtoseconds
	int steps;
	float temperature;
	float cutoff;
	bool thermostat;
	int unit_cells_x;
	int unit_cells_y;
	int unit_cells_z;



	//Material
	float mass;
	float sigma;
	float epsilon;
	float lattice_constant;
	std::string crystal_structure; //bcc,fcc,hcp

	Cell_list* cell_list; 

	//Methods
	void create_cell_list();
	int calculate_number_of_atoms();
	

public:
	void create_list_of_atoms(); //Create all atoms in this class? Convert from fcc to atom positions.
	bool check_input();
	void next_time_step(int current_time_step); //Alter everything in the simulation to get to the next time step.
	void regulate_thermostat(); //Regulate the kinetic energy so that the temperature remains "constant"
	void run_simulation(); //Loop through next_time_step
	void save(); //Save ??? to a .txt file with some structure.
	Simulation (int unit_cells_x,
		int unit_cells_y, // unit_cells is a material parameter.
		int unit_cells_z,
		float time_step,
		int steps,
		float temperature,
		float cutoff,
		float mass,
		float sigma,
		float epsilon,
		float lattice_constant,
		std::string crystal_structure,
		bool thermostat);
	~Simulation ();
	void update_atoms(); // Run through list_of_atoms and .update_atom
	void scc_structure();
	void scc_structure_x(int,int);
	void fcc_structure();
	void fcc_structure_x(int,int);
	void bcc_structure();
	void bcc_structure_x(int,int);
	
	// abort simulation button?
	// handle end of simulation?
};
#endif
