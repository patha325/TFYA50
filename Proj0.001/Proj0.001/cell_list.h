#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <vector>
#include <map>
#include "cell.h"

class Cell_list{
private:
	//Variables
	float cutoff; //Potential cutoff
	float cell_length_x; //Length of side of cell
    float cell_length_y; //Length of side of cell
    float cell_length_z; //Length of side of cell
    float lattice_constant;
    float bulk_length_x;
    float bulk_length_y;
    float bulk_length_z;
	bool pbc_z; //True if periodic boundary condition IS used (I.e. bulk simulation)
	std::vector<Cell*> list_of_cells; //Vector containing all cells
    std::map<int,std::vector<Cell*>> number_to_cell_vector_map; //Vector with neighbouring cells
    
	//Methods
	
    void add_to_map(int,std::vector<Cell*>);

public:
	//Methods
	Cell_list(float cutoff, int unit_cells_x, int unit_cells_y, int unit_cells_z, float lattice_constant, bool new_pbc_z); // unit_cells is a material parameter.
	~Cell_list();
	void add_atom_to_cells(Atom*);
	//int my_sign(float);
    void create_cells();
	std::vector<Atom*> get_neighbours(Atom* atom); 
	void clear_cells();
};

#endif