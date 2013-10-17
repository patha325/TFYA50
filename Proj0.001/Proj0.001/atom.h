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
	Vec next_velocity;
	int cell_number;
	int number_of_neighbours;
		

public:
	//Constructor
	Atom (Vec);
	~Atom ();

	//Getters
	Vec get_velocity();
	Vec get_position();
	int get_cell_number();

	//Setters
	void set_velocity(Vec);
	void set_position(Vec);
	void set_cell_number(int);

	//Other functions
	Vec calculate_force(std::vector<Atom*>);
	float calculate_potential(std::vector<Atom*>);
	Vec distance_vector(Atom*); // Take care of periodic boundry conditions? 
	void next_time_step(); //Alter everything in the atom to get to the next time step.
	void update_atom();	
};

#endif