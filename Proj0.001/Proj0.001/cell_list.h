#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <vector>
#include "cell.h"

class Cell_list{
private:
	//Variables
	float cutoff; //Potential cutoff
	float cell_length; //Length of side of cell
	std::vector<Cell*> list_of_cells; //Vector containing all cells

	//Methods
	void clear_cells(); 
	void create_cell();

public:
	//Methods
	Cell_list(float cutoff); 
	~Cell_list();
	void add_atoms_to_cells(std::vector<Atom*>);
	std::vector<Atom*> get_neighbours(Atom* atom); 
};

#endif