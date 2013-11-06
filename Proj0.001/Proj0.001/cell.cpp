#include "cell.h"

using namespace std;

/*-----------------
CONSTRUCTOR
Sets parameters cell_number
and origin.
-----------------*/
Cell::Cell(int new_cell_number, Vec new_origin){

    cell_number = new_cell_number;
    origin_of_cell = new_origin;
	number_of_atoms_in_cell = 0;
}

/*----------------
DESTRUCTOR
Does nothing
----------------*/
Cell::~Cell(){

    clear_cell();
}


/*---------------
 GETTERS
 --------------*/
int Cell::get_cell_number(){

    return cell_number;
}

Vec Cell::get_origin_of_cell(){

    return origin_of_cell;
}

int Cell::get_number_of_atoms_in_cell(){

    return number_of_atoms_in_cell;
}

std::vector<Atom*> Cell::get_atoms_in_cell(){

    return atoms_in_cell;
}

/*------------------
 SETTERS
-------------------*/
void Cell::set_cell_number(int new_cell_number){

    cell_number = new_cell_number;
}
void Cell::set_origin_of_cell(Vec new_origin){

    origin_of_cell = new_origin;
}

/*--------------------
FUNCTION: add_atom(Atom*)
PARAMETERS: Atom*
RETURNS: Nothing
-
Adds pointer to an Atom 
to the cell.
---------------------*/
void Cell::add_atom(Atom* atom){

    atoms_in_cell.insert(atoms_in_cell.end(), atom);
    number_of_atoms_in_cell++;
}


/*---------------------
FUNCTION: clear_cell()
PARAMETERS: None
RETURNS: Nothing
-
Clears the cell of all 
pointers to Atoms.
---------------------*/
void Cell::clear_cell(){

    atoms_in_cell.clear();
    number_of_atoms_in_cell=0;
}




























