#include "atom.h"

using namespace std;

/*---------------------
CONSTRUCTOR
Parameters: Vec starting_position
Sets starting position
----------------------*/
Atom::Atom(Vec starting_position){

	position = starting_position;
}

/*--------------------
DESTRUCTOR
Does nothing
---------------------*/
Atom::~Atom(){}

/*----------------------
FUNCTION: calculate_force
Paramteters: vector<Atom*>
Returns: Vec (force vector)
-
Calculates force on the atom
from all neighbouring atoms 
within cutoff.
----------------------*/
Vec Atom::calculate_force(vector<Atom*> neigbouring_atoms){

	Vec tmp(0,0,0);
	return tmp;
}
	

/*----------------------
FUNCTION: calculate_potential
Parameters: vector<Atom*>
Returns: float (scalar potential)
-
Calculates (LJ) potential on 
the atom, from all neighbouring
atoms within cutoff.
----------------------*/
float Atom::calculate_potential(vector<Atom*> neighbouring_atoms){

	return 0.5;
}
	
/*----------------------
FUNCTION:distance_vector
Parameters: Atom*
Returns: Vec (vector between atoms)
-
Returns a Vec (vector) which is the 
cartesian vector between the atom
and the parameter atom.
----------------------*/
Vec Atom::distance_vector(Atom* other_atom){

	// Take care of periodic boundry conditions
	Vec tmp(0,0,0);
	return tmp;
}

/*-------------------------
FUNCTION: next_time_step()
Parameters: None
Returns: Nothing
-
Calculates the atom paramters
for the next time step
--------------------------*/
void Atom::next_time_step(){} 

/*--------------------------
FUNCTION: update_atom
Parameters: None
Returns: Nothing
-
Changes the atom paramters
from the state at time t to
time t+1.
--------------------------*/
void Atom::update_atom(){}


// ------- GETTERS --------
Vec Atom::get_velocity(){

	Vec tmp(0,0,0);
	return tmp;
}

Vec Atom::get_position(){

	Vec tmp(0,0,0);
	return tmp;
}

int Atom::get_cell_number(){

	return cell_number;
}


// -------- SETTERS --------
void Atom::set_velocity(Vec newVelocity){}
	
void Atom::set_position(Vec newPosition){}

void Atom::set_cell_number(int){}




