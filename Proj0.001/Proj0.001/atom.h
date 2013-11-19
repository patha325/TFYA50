#ifndef ATOM_H 
#define ATOM_H

#include<vector>
#include <map>

#include "Vec.h"

class Atom{
private:
	Vec position;
	Vec velocity;
	Vec prev_position;
	Vec next_position;
	Vec acceleration;			// Needed to calculate nex pos.
	Vec prev_acceleration;		// Needed to calculate nex pos.
	Vec next_acceleration;		// Needed to calculate nex vel.
	int cell_number;
	int number_of_neighbours;
	float cutoff;
	float bulk_length_x;
	float bulk_length_y;
	float bulk_length_z;
	float lattice_constant;
	float sigma;
	float epsilon;
	float mass;
	float time_step;
	float total_energy;
	float initial_velocity_modulus;
	bool pbc_z;

	Vec distance_vector_pbc(Atom*);
	Vec distance_vector_no_pbc(Atom*);
		
public:
	//Constructor


	Atom (Vec starting_position, Vec new_prev_acceleration, float start_cutoff, int unit_cells_x, int unit_cells_y, int unit_cells_z, float new_lattice_constant,
	float new_sigma, float new_epsilon, float new_mass, float new_time_step, float initial_velocity_modulus, bool new_pbc_z);


	~Atom ();

	//Getters
	Vec get_velocity();
	Vec get_position();
	Vec get_next_position();
	Vec get_acceleration();
	int get_cell_number();
	Vec get_prev_position();

	//Setters
	void set_velocity(Vec);
	void set_position(Vec);
	void set_next_position(Vec);
	void set_prev_position(Vec);
	void set_acceleration(Vec);
	void set_prev_acceleration(Vec);
	void set_cell_number(int);
	void set_cutoff(float);

	//Other functions
	Vec calculate_force(std::vector<Atom*>);
	Vec calculate_acceleration(std::vector<Atom*>);
	float calculate_potential(std::vector<Atom*>);
	Vec calculate_velocity();
	float calculate_kinetic_energy();
	float calculate_temperature(float);
	Vec distance_vector(Atom*); // Take care of periodic boundry conditions? done
	Vec calculate_next_position();
	Vec generate_random_vector();
	void update_atom();	
	vector<Atom*> reduce_neighbours_list(vector<Atom*> original_list);
};

#endif