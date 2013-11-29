#include "atom.h"
#include "vec.h"
#include <math.h>
#include <map>
#include <string>

using namespace std;

/*---------------------
CONSTRUCTOR
Parameters: Vec starting_position
Sets starting position
----------------------*/


Atom::Atom(Vec starting_position, float start_cutoff, int unit_cells_x, int unit_cells_y, int unit_cells_z, float new_lattice_constant,
	float new_sigma, float new_epsilon, float new_mass, float new_time_step, float initial_velocity_modulus, bool new_pbc_z){
	

	position = starting_position;
	prev_acceleration = Vec (0,0,0);
	velocity = initial_velocity_modulus * generate_random_vector(); //generates normed vector in some random direction
	prev_velocity = velocity;
	prev_position = position;
	acceleration = Vec (0,0,0);
	initial_velocity = velocity;
	pbc_z = new_pbc_z;

	cutoff = start_cutoff;
	lattice_constant = new_lattice_constant;
	sigma = new_sigma;
	epsilon = new_epsilon;
	bulk_length_x = unit_cells_x*lattice_constant;
	bulk_length_y = unit_cells_y*lattice_constant;
	bulk_length_z = unit_cells_z*lattice_constant;
	mass = new_mass;
	time_step = new_time_step;
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
Vec Atom::calculate_force(Vec distance, float distance_length){

	Vec tmp_force (0,0,0);
	float q = sigma/distance_length;
	tmp_force =(48/distance_length*epsilon*(pow(q, 12)-0.5f*pow(q, 6)))*distance.normalize();

	return tmp_force;

}

/*----------------------
FUNCTION: calculate_pressure
Paramteters: vector<Atom*>
Returns: float (pressure)
-
Calculates pressure on the atom
from all neighbouring atoms 
within cutoff.
----------------------*/
float Atom::calculate_pressure(Vec tmp_force, float distance_length){

    float tmp_pressure = 0;

	tmp_pressure = tmp_force.length()*distance_length;

	return tmp_pressure;
}

/*----------------------
FUNCTION: calculate__and_set_acceleration
Paramteters: vector<Atom*>
Returns: Vec (acceleration vector)
-
Calculates and sets acceleration on the atom
from closest atom.
----------------------*/

void Atom::calculate_and_set_acceleration(Vec force){

	acceleration = (-1/mass)*force;
	
	return;
}

/*----------------------
FUNCTION: calculate_potential
Parameters: vector<Atom*>
Returns: float (scalar potential)
-
Calculates (LJ) potential on 
the atom, from closest atom.
----------------------*/
float Atom::calculate_potential(float distance_length, Atom* other_atom){

	float tmp_potential = 0;
	float q = sigma/distance_length; 
	tmp_potential = 4*epsilon*(pow(q,12)-pow(q,6));

	/*
	if(atom_number == 0){
		cout << "Distance to atom #" << other_atom->get_atom_number() << " is: " <<distance_length << endl;
	}
	*/

	return tmp_potential;
}

/*----------------------
FUNCTION: calculate_and_set_velocity
Parameters: none
Returns: Vec (velocity)
-
Calculates and sets velocity on 
the atom.
----------------------*/

void Atom::calculate_and_set_velocity(){

	velocity.setCoords(
		prev_velocity.getX() + (time_step/2)*(acceleration.getX() + prev_acceleration.getX()),
		prev_velocity.getY() + (time_step/2)*(acceleration.getY() + prev_acceleration.getY()),
		prev_velocity.getZ() + (time_step/2)*(acceleration.getZ() + prev_acceleration.getZ()));

	return;
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

	float k_b = 8.617342e-5f; //[eV][K]^{-1}
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
	
	if(pbc_z) return distance_vector_pbc(other_atom);
	else return distance_vector_no_pbc(other_atom);
	
/*
	Vec v1 = distance_vector_pbc(other_atom);
	Vec v2 = distance_vector_no_pbc(other_atom);

	if(position.getZ()>22.0){
		cout << "----------------------" << endl;
		cout << "Reference atom position: " << position << endl;
		cout << "Other atom position: " << other_atom->get_position() << endl;
		cout << "Same bulk:" << other_atom->get_position() - position << endl;
		cout << "V1: " << v1 << endl;
		cout << "V2: " << v2 << endl;
		system("pause");
		cout << "----------------------" << endl << endl;
	}
	if (v1!=v2) cout << "ERROR: Does not give the same distance vector" << endl;
*/

}

/*--------------------------------------------
Function for distance vector when using PBC in
all directions.
--------------------------------------------*/
Vec Atom::distance_vector_pbc(Atom* other_atom){
	Vec tmp = other_atom->position;

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
	/*if(shortest_vec.length() < 0.5 && atom_number != other_atom->get_atom_number()){
		cout << "Atom " << atom_number << " mkt nära" << other_atom->get_atom_number() << endl;
	}*/
	// Returns the vector to the closest atom from list
	return shortest_vec;
}

/*----------------------------------------
Function for distance vector when NOT using 
PBC in all directions.
-----------------------------------------*/
Vec Atom::distance_vector_no_pbc(Atom* other_atom){

//	cout << "----------------------" << endl;

	Vec other_atom_vector = other_atom->get_position();
	Vec saved_distance = other_atom_vector - position;
	Vec new_atom_vector;

//	cout << "Reference atom position: " << position << endl;


	for(int x = -1; x <= 1; x++){
		for(int y = -1; y <= 1; y++){
			for(int z = -1; z <= 0; z++){
				new_atom_vector = other_atom_vector+Vec(bulk_length_x*x, bulk_length_y*y, bulk_length_z*z);
				Vec tmp_distance = new_atom_vector - position;

/*				
				cout << "Other atom position: " << new_atom_vector << endl;
				cout << "Distance vector: " << tmp_distance << endl;
				cout << "Distance: " << tmp_distance.length() << endl;
*/
				if(tmp_distance.length()<saved_distance.length()){
//					cout << "- Saved -" << endl;
					saved_distance = tmp_distance;
				}
			}
		}
	}

/*
	cout << "FINAL DISTANCE: " << saved_distance << endl;
	cout << "---------------------------" << endl << endl;
	system("pause");
*/

	return saved_distance;
}

/*-------------------------
FUNCTION: calculate_and_set_position()
Parameters: None
Returns: None
-
Calculates and sets position for an atom
--------------------------*/
void Atom::calculate_and_set_position(){
	float time_step2 = time_step*time_step;
	//You should not change atom attributes, Markus!
	Vec next_position (0,0,0);

	next_position.setCoords(
		prev_position.getX() + prev_velocity.getX()*time_step + prev_acceleration.getX()*1/2*time_step2,
		prev_position.getY() + prev_velocity.getY()*time_step + prev_acceleration.getY()*1/2*time_step2, 
		prev_position.getZ() + prev_velocity.getZ()*time_step + prev_acceleration.getZ()*1/2*time_step2);


	//Check if atom outside cells and move it inside if it is. 
	//If surface simulation, no atoms are moved like this.
	if(pbc_z){

		//x
		while(next_position.getX() < 0 || next_position.getX() >= bulk_length_x){
			int sign = my_sign(next_position.getX());
			next_position.setX(next_position.getX()-sign*bulk_length_x);
		}
		//y
		while(next_position.getY() < 0 || next_position.getY() >= bulk_length_y){
			int sign = my_sign(next_position.getY());
			next_position.setY(next_position.getY()-sign*bulk_length_y);
		}
		//z
		while(next_position.getZ() < 0 || next_position.getZ() >= bulk_length_z){
			int sign = my_sign(next_position.getZ());
			next_position.setZ(next_position.getZ()-sign*bulk_length_z);
		}
	}

	/*position = next_position;
	if((abs(position.getX()-prev_position.getX()) > 0.5 && abs(position.getX()-prev_position.getX()) < 25.5)||
		(abs(position.getY()-prev_position.getY()) > 0.5 && abs(position.getY()-prev_position.getY()) < 25.5)||
		(abs(position.getZ()-prev_position.getZ()) > 0.5 && abs(position.getZ()-prev_position.getZ()) < 25.5))
	{
			cout << "hopp i position, atom " << atom_number << endl;
			cout << "  position:     " << position << endl;
			cout << "  prev_position " << prev_position << endl << endl;

	}*/

	position = next_position;
	
	return;
	// Change the cell number? Should there be a call for that? Has been added in add_atoms_to_cell in cell_list
}

/* ------------------------------
FUNCTION: sign
PARAMETERS: float
RETURN: int
-
Returns +1 or -1 depending on the
sign of the incoming float
------------------------------*/
int Atom::my_sign(float number){
	int sign;
	if(number < 0){
		sign = -1;
	}
	else{
		sign = +1;
	}
	return sign;
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
	float x = ((float) rand() / (RAND_MAX));
	float y = ((float) rand() / (RAND_MAX));
	float z = ((float) rand() / (RAND_MAX));
	float vec_mod = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	return Vec (x/vec_mod, y/vec_mod, z/vec_mod);
}

/*--------------------------
FUNCTION: update_atom
Parameters: None
Returns: Nothing
-
Changes position of atom
--------------------------*/
void Atom::update_atom(){

	prev_position = position;
	prev_velocity = velocity;
	prev_acceleration = acceleration;
}



void Atom::update_neighbour_list(vector<Atom*> new_neighbours){

	atom_neighbours = new_neighbours;
	number_of_neighbours = new_neighbours.size();
	//cout << "update_neighbour_list atom_neighbours.size(): " << atom_neighbours.size() << endl;
}

void Atom::clear_tmp_force(){

	tmp_force = Vec (0,0,0);
	return;
}

// ------- GETTERS --------
int Atom::get_atom_number(){ return atom_number; }
Vec Atom::get_velocity(){ return velocity; }
Vec Atom::get_position(){ return position; }
Vec Atom::get_acceleration(){ return acceleration; }
Vec Atom::get_prev_acceleration(){ return prev_acceleration; }
int Atom::get_cell_number(){ return cell_number; }
Vec Atom::get_prev_position(){ return prev_position; }
Vec Atom::get_initial_velocity() { return initial_velocity; }
Vec Atom::get_initial_position() { return initial_position; }
float Atom::get_mass(){ return mass; }
vector<Atom*> Atom::get_atom_neighbours(){ return atom_neighbours; }
Vec Atom::get_tmp_force(){ return tmp_force; }
int Atom::get_number_of_neighbours(){return number_of_neighbours;}
Vec Atom::get_prev_velocity(){return prev_velocity;}



// -------- SETTERS --------
void Atom::set_atom_number(int new_atom_number){ atom_number = new_atom_number; }
void Atom::set_velocity(Vec new_velocity){ velocity = new_velocity; }
void Atom::set_cell_number(int new_cell_number){ 
	//if(new_cell_number!=cell_number) cout << "Moved atom " << atom_number << endl;
	cell_number = new_cell_number; 
}
void Atom::set_initial_velocity(Vec new_initial_velocity) { initial_velocity = new_initial_velocity; }
void Atom::set_initial_position(Vec new_initial_position) { initial_position = new_initial_position; }
void Atom::set_cutoff(float new_cutoff){ cutoff = new_cutoff; }
void Atom::add_tmp_force(Vec new_tmp_force){ tmp_force = tmp_force + new_tmp_force; }







