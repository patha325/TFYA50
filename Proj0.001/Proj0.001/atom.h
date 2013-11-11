#ifndef ATOM_H 
#define ATOM_H

#include<vector>

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
	float sigma;
	float epsilon;
	float mass;
	float time_step;
	float total_energy;
	float initial_velocity_modulus;
		

public:
	//Constructor
	Atom (Vec,Vec,float,float,float,float,float,float,float,float,float);
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
};

#endif