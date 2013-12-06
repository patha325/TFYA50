#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <map>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <time.h>

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
	float pressure;
	float cutoff;
	bool thermostat;
	bool equilibrium;
	bool pbc_z;
	int thermostat_update_freq;
	bool old_sim;
	float initial_velocity_modulus;
	int unit_cells_x;
	int unit_cells_y;
	int unit_cells_z;
	float volume;
	float total_energy;
	float Diff_coeff;
	float self_diff_coeff;
	float last_MSD;
	float prev_diff_coeff;
	map<int, vector<Vec>> atom_positions;
	int eq_time_steps;
	//Boltzmann constant
	float k_b;
	//Planck's constant
	float hbar;
	bool save_atom_positions;



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
	//Methods
	void create_list_of_atoms(); //Create all atoms in this class? Convert from fcc to atom positions.
	bool check_input();
	void next_time_step(int current_time_step); //Alter everything in the simulation to get to the next time step.
	void calculate_and_set_velocity(Atom* atom,double current_time_step); //Sets velocity. Checks if system has thermostat or not
	float calculate_specific_heat();
	float calculate_MSD(Atom*);
	Vec distance_between_now_and_initial(Atom*);
	void regulate_thermostat(); //Regulate the kinetic energy so that the temperature remains "constant"
	void run_simulation(); //Loop through next_time_step and return last state
	void save(); //Save ??? to a .txt file with some structure.
	void update_atoms(); // Run through list_of_atoms and .update_atom
	void scc_structure();
	void scc_structure_x(int,int);
	void scc_corrector();
	void fcc_structure();
	void fcc_structure_x(int,int);
	void fcc_corrector();
	void bcc_structure();
	void bcc_structure_x(int,int);
	void bcc_corrector();
	std::vector<Atom*> get_list_of_atoms();
	int get_number_of_atoms();
	void end_of_simulation();
	void read_old_sim();

	//Constructors
	//std::ofstream fs2;
	Simulation (
		int unit_cells_x,
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
		bool thermostat,
		bool equilibrium,
		bool pbc_z,
		int thermostat_update_freq,
		bool old_sim,
		bool save_atom_positions);
	Simulation(Simulation* old_simulation, int steps, bool equilibrium); //Take off where we left off-constructor

	//Destructor
	~Simulation();

	//Getters
	float get_time_step();
	int get_steps();
	float get_temperature();
	float get_cutoff();
	bool get_thermostat();
	bool get_pbc_z();
	float get_initial_velocity_modulus();
	int get_unit_cells_x();
	int get_unit_cells_y();
	int get_unit_cells_z();
	float get_total_energy();
	float get_mass();
	float get_sigma();
	float get_epsilon();
	float get_lattice_constant();
	std::string get_crystal_structure();
	Cell_list* get_cell_list(); 
	bool get_save_atom_positions();

	void configure_data(int);
	
	// abort simulation button?
	// handle end of simulation?
};
#endif
