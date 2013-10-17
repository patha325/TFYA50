#ifndef ATOM_H 
#define ATOM_H

#include<vector>

class Atom{
private:

	std::vector<float> position;
	std::vector<float> velocity;
	std::vector<float> prev_position;
	std::vector<float> next_position;
	std::vector<float> next_velocity;
	int cell_number;
	int number_of_neighbours;
		

public:
	std::vector<float> calculate_force(std::vector<Atom*>);
	float calculate_potential(std::vector<Atom*>);
	std::vector<float> distance(Atom*); // Take care of periodic boundry conditions? 
	std::vector<float> get_velocity();
	void set_velocity(std::vector<float>);
	std::vector<float> get_position();
	void set_position(std::vector<float>);
	int get_cell_number();
	void set_cell_number(int);
	void next_time_step(); //Alter everything in the atom to get to the next time step.
	Atom (std::vector<float> position);
	~Atom ();
	void update_atom();




};

#endif