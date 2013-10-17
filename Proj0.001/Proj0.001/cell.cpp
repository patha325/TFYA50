#include "cell.h"

using namespace std;

/*
int cell_number; //Do we need this?
std::vector<float> origin_of_cell;// Lower left corner. ALWAYS
std::vector<Atom*> atoms_in_cell;
int number_of_atoms_in_cell;
*/

/*-----------------
CONSTRUCTOR
Sets parameters cell_number
and origin.
-----------------*/
Cell::Cell(int cell_number, std::vector<float> origin){}

/*----------------
DESTRUCTOR
Does nothing
----------------*/
Cell::~Cell(){}


/*--------------------
FUNCTION: add_atom(Atom*)
PARAMETERS: Atom*
RETURNS: Nothing
-
Adds pointer to an Atom 
to the cell.
---------------------*/
void Cell::add_atom(Atom* atom){}


/*---------------------
FUNCTION: clear_cell()
PARAMETERS: None
RETURNS: Nothing
-
Clears the cell of all 
pointers to Atoms.
---------------------*/
void Cell::clear_cell(){}
