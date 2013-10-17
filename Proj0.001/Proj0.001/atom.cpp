#include "atom.h"

using namespace std;

vector<float> Atom::calculate_force(vector<Atom*> neigbouring_atoms){

	vector<float> tmp;
	tmp.insert(tmp.begin(), 0.5);
	tmp.insert(tmp.begin(), 0.6);
	tmp.insert(tmp.begin(), 0.7);
	return tmp;
}
	

float Atom::calculate_potential(vector<Atom*> neighbouring_atoms){

	return 0.5;
}
	
vector<float> Atom::distance(Atom*){

	// Take care of periodic boundry conditions
	vector<float> tmp;
	tmp.insert(tmp.begin(), 0.5);
	tmp.insert(tmp.begin(), 0.6);
	tmp.insert(tmp.begin(), 0.7);
	return tmp;
} 

vector<float> Atom::get_velocity(){

	vector<float> tmp;
	tmp.insert(tmp.begin(), 0.5);
	tmp.insert(tmp.begin(), 0.6);
	tmp.insert(tmp.begin(), 0.7);
	return tmp;
}

void Atom::set_velocity(std::vector<float>){}
	
vector<float> Atom::get_position(){

	vector<float> tmp;
	tmp.insert(tmp.begin(), 0.5);
	tmp.insert(tmp.begin(), 0.6);
	tmp.insert(tmp.begin(), 0.7);
	return tmp;
}

void Atom::set_position(vector<float>){

}

int Atom::get_cell_number(){

	return 0;
}

void Atom::set_cell_number(int){}

void Atom::next_time_step(){} //Alter everything in the atom to get to the next time step.

Atom::Atom(vector<float> starting_position){

	position = starting_position;
}

Atom::~Atom(){}

void Atom::update_atom(){}
