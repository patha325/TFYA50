#include <iostream>
#include "cell_list.h"
#include "cell.h"
#include "atom.h"

using namespace std;

/* ------------------------------
CONSTRUCTOR
PARAMETERS: float cutoff [cutoff distance]
-
Creates the list_of_cells vector
and the cells in the vector.
------------------------------ */
Cell_list::Cell_list(float cutoff){} 


/* ------------------------------
DESTRUCTOR
PARAMETERS: None
-
Destroys cells in list_of_cells
------------------------------ */
Cell_list::~Cell_list(){}


/* ---
PUBLIC
--- */

/* ------------------------------
FUNCTION: Cell_list::add_atoms_to_cells()
PARAMETERS: vector<Atoms*>
RETURN: void
-
Adds all Atoms to the right cells
------------------------------ */
void Cell_list::add_atoms_to_cells(vector<Atom*>){}



/* ------------------------------
FUNCTION: Cell_list::get_neighbours()
PARAMETERS: Atom* atom
RETURN: vector<Atom*>
-
Returns a vector with all Atoms which
are neighbours to parameter Atom.
------------------------------ */
vector<Atom*> Cell_list::get_neighbours(Atom* atom){

	vector<Atom*> tmp;
	return tmp;
} 


/*---
PRIVATE
--- */

/* ------------------------------
FUNCTION: Cell_list::clear_cells()
PARAMETERS: None
RETURN: void
-
Removes all atoms from all cells
------------------------------ */
void Cell_list::clear_cells(){}


/* ------------------------------
FUNCTION: Cell_list::create_cell()
PARAMETERS: None
RETURN: void
-
Creates Cell-objects and adds them
to the list_of_cells vector.
------------------------------ */
void Cell_list::create_cell(){}


