#include <iostream>
#include "cell_list.h"
#include "cell.h"
#include "atom.h"

using namespace std;

//private
void Cell_list::clear_cells(){}
void Cell_list::create_cell(){}


//public
Cell_list::Cell_list(float cutoff){} 

Cell_list::~Cell_list(){}

void Cell_list::add_atoms_to_cells(std::vector<Atom*>){}

vector<Atom*> Cell_list::getNeighbours(Atom* atom){

	vector<Atom*> tmp;
	return tmp;
} //Returns vector of atoms in nearest cells, with respect to paramteter atom.
