#include "atom.h"
#include "vec.h"
#include <math.h>

using namespace std;

/*---------------------
CONSTRUCTOR
Parameters: Vec starting_position
Sets starting position
----------------------*/

Atom::Atom(Vec starting_position, Vec new_prev_acceleration, float start_cutoff, float unit_cells_x, float unit_cells_y, float unit_cells_z, 
	float new_sigma, float new_epsilon, float new_mass, float new_time_step, float initial_velocity_modulus){


	position = starting_position;
	prev_acceleration = new_prev_acceleration;
	cutoff = start_cutoff;
	sigma = new_sigma;
	epsilon = new_epsilon;
	bulk_length_x = unit_cells_x*sigma;
	bulk_length_y = unit_cells_y*sigma;
	bulk_length_z = unit_cells_z*sigma;
	mass = new_mass;
	time_step = new_time_step;
	velocity = initial_velocity_modulus * generate_random_vector();
	cout << "velocity " << velocity << endl;
	prev_position = position;
	next_position = Vec (0, 0, 0);
	acceleration = Vec (0, 0, 0);
	next_acceleration = Vec (0, 0, 0);
}

/*--------------------
DESTRUCTOR
Does nothing
---------------------*/
Atom::~Atom(){}

//----------Other functions-------------------


/*----------------------
FUNCTION: calculate_force
Paramteters: vector<Atom*>
Returns: Vec (force vector)
-
Calculates force on the atom
from all neighbouring atoms 
within cutoff.
----------------------*/
Vec Atom::calculate_force(vector<Atom*> neighbouring_atoms){ 
    
	Vec tmp_force (0,0,0);
	for(string::size_type i=0; i < neighbouring_atoms.size(); i++){
		// string::size_type ist för int eftersom .size() returnerar en unsigned int, blir varning annars.

		float r= distance_vector(neighbouring_atoms[i]).length();
		float r2 = pow(r,-13);
		float r3 = pow(r,-7);
		tmp_force += -48*(r2-r3/2)*distance_vector(neighbouring_atoms[i]).normalize();
	}
	return tmp_force;
}

/*----------------------
FUNCTION: calculate_acceleration
Paramteters: vector<Atom*>
Returns: Vec (acceleration vector)
-
Calculates acceleration on the atom
from closest atom.
----------------------*/

Vec Atom::calculate_acceleration(vector<Atom*> neighbouring_atoms){

	Vec tmp_acc (0,0,0);
	Vec tmp_force = calculate_force(neighbouring_atoms);
	tmp_acc.setCoords(-tmp_force.getX()/mass,-tmp_force.getY()/mass,-tmp_force.getZ()/mass);
	
	/*prev_acceleration = acceleration;
	acceleration = tmp_acc; Kanske int här...*/
	
	//cout << "acc " << tmp_acc << endl; 
	return tmp_acc;
}

/*----------------------
FUNCTION: calculate_potential
Parameters: vector<Atom*>
Returns: float (scalar potential)
-
Calculates (LJ) potential on 
the atom, from closest atom.
----------------------*/
float Atom::calculate_potential(vector<Atom*> neighbouring_atoms){

	float tmp_potential = 0;
	for(string::size_type i = 0; i < neighbouring_atoms.size(); i++){
		// string::size_type ist för int eftersom .size() returnerar en unsigned int, blir varning annars.

		Vec closest_vector_tmp = distance_vector(neighbouring_atoms[i]);
		float tmp_distance = closest_vector_tmp.length();
		tmp_potential += 4*epsilon*(pow(sigma/tmp_distance,12)-pow(sigma/tmp_distance,6));
	}
	return tmp_potential;
}

/*----------------------
FUNCTION: calculate_velocity
Parameters: none
Returns: Vec (velocity)
-
Calculates velocity on 
the atom.
----------------------*/

Vec Atom::calculate_velocity(){

	Vec position_diff = position - prev_position;
	Vec tmp_velocity;
	tmp_velocity.setCoords(position_diff.getX()/time_step,position_diff.getY()/time_step,position_diff.getZ()/time_step);
	return tmp_velocity;
}

/*----------------------
FUNCTION: calculate_kinetic_energy
Parameters: none
Returns: float (Kinetic_energy)
-
Calculates kinetik energy on 
the atom.
----------------------*/

float Atom::calculate_kinetic_energy(){
	
	float tmp_kinetic_energy = mass*pow(velocity.length(),2)/2;
	return tmp_kinetic_energy;
}

/*----------------------
FUNCTION: calculate_temperature
Parameters: none
Returns: float (temperature)
-
Calculates temperature on 
the atom. Ganska vagt kanske...
----------------------*/

float Atom::calculate_temperature(float E_kin){
	

	float k_b = 8.617342e-5; //[eV][K]^{-1}
	return (2*E_kin)/(3*k_b);

}
	
/*----------------------
FUNCTION:distance_vector
Parameters: Atom*
Returns: Vec (vector between atoms)
-
Returns a Vec (vector) which is the 
cartesian vector between the atom
and the parameter atom. Jobbar förmodligen i nm när den används i potential mm...
----------------------*/
Vec Atom::distance_vector(Atom* other_atom){
	Vec tmp =other_atom->position;

	/* Check the x, y and z coordinate for both atoms */
	float x1 = position.getX();
	float x2 = tmp.getX();
	float y1 = position.getY();
	float y2 = tmp.getY();
	float z1 = position.getZ();
	float z2 = tmp.getZ();

	/* Check distance inside this periodic bulk structure*/
	Vec l1 = tmp - position;

	/* Check distance to atom in neighbouring periodic bulk structures, only checking the 7 closest */

	Vec l2;
	Vec l3;
	Vec l4;
	Vec l5;
	Vec l6;
	Vec l7;
	Vec l8;
	
	if(x1>=x2 && y1>=y2 && z1>=z2){
		/* Check x + bulk_length, y + bulk_length, z + bulk_length and all combinations */
		l2 = (tmp + Vec(bulk_length_x,0,0)) - position;
		l3 = (tmp + Vec(bulk_length_x,bulk_length_y,0)) - position;
		l4 = (tmp + Vec(bulk_length_x,0,bulk_length_z)) - position;
		l5 = (tmp + Vec(bulk_length_x,bulk_length_y,bulk_length_z)) - position;
		l6 = (tmp + Vec(0,bulk_length_y,0)) - position;
		l7 = (tmp + Vec(0,bulk_length_y,bulk_length_z)) - position;
		l8 = (tmp + Vec(0,0,bulk_length_z)) - position;
	}
	else if(x1>=x2 && y1>=y2 && z1<=z2){
		/* Check x + bulk_length, y + bulk_length, z - bulk_length and all combinations */
		l2 = (tmp + Vec(bulk_length_x,0,0)) - position;
		l3 = (tmp + Vec(bulk_length_x,bulk_length_y,0)) - position;
		l4 = (tmp + Vec(bulk_length_x,0,-bulk_length_z)) - position;
		l5 = (tmp + Vec(bulk_length_x,bulk_length_y,-bulk_length_z)) - position;
		l6 = (tmp + Vec(0,bulk_length_y,0)) - position;
		l7 = (tmp + Vec(0,bulk_length_y,-bulk_length_z)) - position;
		l8 = (tmp + Vec(0,0,-bulk_length_z)) - position;
	}
	else if(x1>=x2 && y1<=y2 && z1>=z2){
		/* Check x + bulk_length, y - bulk_length, z + bulk_length and all combinations */
		l2 = (tmp + Vec(bulk_length_x,0,0)) - position;
		l3 = (tmp + Vec(bulk_length_x,-bulk_length_y,0)) - position;
		l4 = (tmp + Vec(bulk_length_x,0,bulk_length_z)) - position;
		l5 = (tmp + Vec(bulk_length_x,-bulk_length_y,bulk_length_z)) - position;
		l6 = (tmp + Vec(0,-bulk_length_y,0)) - position;
		l7 = (tmp + Vec(0,-bulk_length_y,bulk_length_z)) - position;
		l8 = (tmp + Vec(0,0,bulk_length_z)) - position;

	}
	else if(x1>=x2 && y1<=y2 && z1<=z2){
		/* Check x + bulk_length, y - bulk_length, z - bulk_length and all combinations */
		l2 = (tmp + Vec(bulk_length_x,0,0)) - position;
		l3 = (tmp + Vec(bulk_length_x,-bulk_length_y,0)) - position;
		l4 = (tmp + Vec(bulk_length_x,0,-bulk_length_z)) - position;
		l5 = (tmp + Vec(bulk_length_x,-bulk_length_y,-bulk_length_z)) - position;
		l6 = (tmp + Vec(0,-bulk_length_y,0)) - position;
		l7 = (tmp + Vec(0,-bulk_length_y,-bulk_length_z)) - position;
		l8 = (tmp + Vec(0,0,-bulk_length_z)) - position;
	}
	else if(x1<=x2 && y1>=y2 && z1>=z2){
		/* Check x - bulk_length, y + bulk_length, z + bulk_length and all combinations */
		l2 = (tmp + Vec(-bulk_length_x,0,0)) - position;
		l3 = (tmp + Vec(-bulk_length_x,bulk_length_y,0)) - position;
		l4 = (tmp + Vec(-bulk_length_x,0,bulk_length_z)) - position;
		l5 = (tmp + Vec(-bulk_length_x,bulk_length_y,bulk_length_z)) - position;
		l6 = (tmp + Vec(0,bulk_length_y,0)) - position;
		l7 = (tmp + Vec(0,bulk_length_y,bulk_length_z)) - position;
		l8 = (tmp + Vec(0,0,bulk_length_z)) - position;
	}
	else if(x1<=x2 && y1<=y2 && z1>=z2){
		/* Check x - bulk_length, y - bulk_length, z + bulk_length and all combinations */
		l2 = (tmp + Vec(-bulk_length_x,0,0)) - position;
		l3 = (tmp + Vec(-bulk_length_x,-bulk_length_y,0)) - position;
		l4 = (tmp + Vec(-bulk_length_x,0,bulk_length_z)) - position;
		l5 = (tmp + Vec(-bulk_length_x,-bulk_length_y,bulk_length_z)) - position;
		l6 = (tmp + Vec(0,-bulk_length_y,0)) - position;
		l7 = (tmp + Vec(0,-bulk_length_y,bulk_length_z)) - position;
		l8 = (tmp + Vec(0,0,bulk_length_z)) - position;
	}
	else if(x1<=x2 && y1>=y2 && z1<=z2){
		/* Check x - bulk_length, y + bulk_length, z - bulk_length and all combinations */
		l2 = (tmp + Vec(-bulk_length_x,0,0)) - position;
		l3 = (tmp + Vec(-bulk_length_x,bulk_length_y,0)) - position;
		l4 = (tmp + Vec(-bulk_length_x,0,-bulk_length_z)) - position;
		l5 = (tmp + Vec(-bulk_length_x,bulk_length_y,-bulk_length_z)) - position;
		l6 = (tmp + Vec(0,bulk_length_y,0)) - position;
		l7 = (tmp + Vec(0,bulk_length_y,-bulk_length_z)) - position;
		l8 = (tmp + Vec(0,0,-bulk_length_z)) - position;
	}
	else{
		/* Check x - bulk_length, y - bulk_length, z - bulk_length and all combinations */
		l2 = (tmp - Vec(bulk_length_x,0,0)) - position;
		l3 = (tmp - Vec(bulk_length_x,bulk_length_y,0)) - position;
		l4 = (tmp - Vec(bulk_length_x,0,bulk_length_z)) - position;
		l5 = (tmp - Vec(bulk_length_x,bulk_length_y,bulk_length_z)) - position;
		l6 = (tmp - Vec(0,bulk_length_y,0)) - position;
		l7 = (tmp - Vec(0,bulk_length_y,bulk_length_z)) - position;
		l8 = (tmp - Vec(0,0,bulk_length_z)) - position;
	}


	/* Check which one of the 8 atoms that are closest*/
	float shortest_distance = l1.length();
	Vec shortest_vec = l1;

	if(l2.length() < shortest_distance){
		shortest_distance = l2.length();
		shortest_vec = l2;
	}
	if(l3.length() < shortest_distance){
		shortest_distance = l3.length();
		shortest_vec = l3;
	}
	if(l4.length() < shortest_distance){
		shortest_distance = l4.length();
		shortest_vec = l4;
	}
	if(l5.length() < shortest_distance){
		shortest_distance = l5.length();
		shortest_vec = l5;
	}
	if(l6.length() < shortest_distance){
		shortest_distance = l6.length();
		shortest_vec = l6;
	}
	if(l7.length() < shortest_distance){
		shortest_distance = l7.length();
		shortest_vec = l7;
	}
	if(l8.length() < shortest_distance){
		shortest_distance = l8.length();
		shortest_vec = l8;
	}

	// Returns the vector to the closest atom from list
	return shortest_vec;
}

/*-------------------------
FUNCTION: calculate_next_position()
Parameters: None
Returns: next_position
-
Calculates next position for an atom
--------------------------*/
Vec Atom::calculate_next_position(){
	float time_step2 = time_step*time_step;
	//You should not change atom attributes, Markus!
	Vec next_position (0,0,0);

	next_position.setCoords(
		position.getX() + velocity.getX()*time_step + acceleration.getX()*1/2*time_step2,
		position.getY() + velocity.getY()*time_step + acceleration.getY()*1/2*time_step2, 
		position.getZ() + velocity.getZ()*time_step + acceleration.getZ()*1/2*time_step2);


	return next_position;
	// Change the cell number? Should there be a call for that? Has been added in add_atoms_to_cell in cell_list
}

/*--------------------------
FUNCTION: generate_random_vector
Parameters: None
Returns: Vec
-
Generates a random Vec but
with modulus 1
--------------------------*/

Vec Atom::generate_random_vector(){
	float x = ((double) rand() / (RAND_MAX));
	float y = ((double) rand() / (RAND_MAX));
	float z = ((double) rand() / (RAND_MAX));
	float vec_mod = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	return Vec (x/vec_mod, y/vec_mod, z/vec_mod);
}

/*--------------------------
FUNCTION: update_atom
Parameters: None
Returns: Nothing
-
Changes the atom paramters
from the state at time t to
time t+time_step.
Acceleration and prev_acceleration 
is set when calculate_acceleration 
is called.
--------------------------*/
void Atom::update_atom(){ // Needed??

}


// ------- GETTERS --------
Vec Atom::get_velocity(){
	
	return velocity;
}

Vec Atom::get_position(){

	return position;
}

Vec Atom::get_next_position(){
	return next_position;
}

Vec Atom::get_acceleration(){

	return acceleration;
}

int Atom::get_cell_number(){

	return cell_number;
}

Vec Atom::get_prev_position(){

	return prev_position;
}


// -------- SETTERS --------
void Atom::set_velocity(Vec new_velocity){
	velocity = new_velocity;
	return;
}
	
void Atom::set_position(Vec new_position){
	position = new_position;
	return;
}

void Atom::set_cell_number(int new_cell_number){
	cell_number = new_cell_number;
	return;
}

void Atom::set_prev_position(Vec new_position){

	prev_position = new_position;
	return;
}

void Atom::set_next_position(Vec new_next_position){

	next_position = new_next_position;
	return;
}

void Atom::set_acceleration(Vec new_acceleration){
	
	acceleration = new_acceleration;
	return;
}

void Atom::set_prev_acceleration(Vec new_acceleration){

	prev_acceleration = new_acceleration;
	return;
}

void Atom::set_cutoff(float new_cutoff){
	cutoff = new_cutoff;
	return;
}




