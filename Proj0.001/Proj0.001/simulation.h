#ifndef SIMULATION_H
#define SIMULATION_H

#include<vector>
#include<string>
#include "atom.h"
#include "cell_list.h"
class Simulation{
private:
	//Variables
	std::vector<Atom*> list_of_atoms;
	int number_of_atoms;
	int time_step; //Given in femtoseconds
	int steps;
	float temperature;
	float cutoff;
	bool thermostat;

	//Material
	float mass;
	float sigma;
	float epsilon;
	float lattice_constant;
	string crystal_structure; //bcc,fcc,hcp

	Cell_list* cell_list; 

	//Methods
	void create_cell_list();
	

public:
	void create_list_of_atoms(); //Create all atoms in this class? Convert from fcc to atom positions.
	void next_time_step(); //Alter everything in the simulation to get to the next time step.
	void regulate_thermostat(); //Regulate the kinetic energy so that the temperature remains "constant"
	void run_simulation(); //Loop through next_time_step
	void save(); //Save ??? to a .txt file with some structure.
	Simulation (int number_of_atoms,
		int time_step,
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

};
#endif