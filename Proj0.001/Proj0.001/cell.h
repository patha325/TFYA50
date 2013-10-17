#ifndef CELL_H
#define CELL_H

#include<vector>
#include"atom.h"

class Cell{
private:
	int cell_number; //Do we need this?
	std::vector<float> origin_of_cell;// Lower left corner. ALWAYS
	std::vector<Atom*> atoms_in_cell;
	int number_of_atoms_in_cell;

public:
	Cell(int cell_number, std::vector<float> origin);
	~Cell();
	void add_atom(Atom*);
	void clear_cell();
};

#endif