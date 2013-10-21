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
	Vec tmp =other_atom->position;
	tmp-=position;
	// Take care of periodic boundry conditions
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
void Atom::update_atom(){
position = next_position;
velocity = next_velocity;
// Change the cell number? Should there be a call for that?
}


// ------- GETTERS --------
Vec Atom::get_velocity(){
	
	return velocity;
}

Vec Atom::get_position(){

	return position;
}

int Atom::get_cell_number(){

	return cell_number;
}


// -------- SETTERS --------
void Atom::set_velocity(Vec newVelocity){}
	
void Atom::set_position(Vec newPosition){}

void Atom::set_cell_number(int){}




