#ifndef CELL_H
#define CELL_H

#include <vector>
#include "atom.h"
#include "vec.h"

class Cell{
private:
	int cell_number; //Do we need this?
	Vec origin_of_cell;// Lower left corner. ALWAYS
	std::vector<Atom*> atoms_in_cell;
	int number_of_atoms_in_cell;

public:
	Cell(int cell_number, Vec);
	~Cell();
	void add_atom(Atom*);
	void clear_cell();
};

#endif