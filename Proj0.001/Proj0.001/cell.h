#ifndef CELL_H
#define CELL_H

#include <vector>
#include <map>
#include "atom.h"
#include "vec.h"

class Cell{
private:
	int cell_number; //Do we need this? Yes we do!
	Vec origin_of_cell;// Lower left corner. ALWAYS
	std::vector<Atom*> atoms_in_cell;
	int number_of_atoms_in_cell;

public:
    //Constructors
	Cell(int cell_number, Vec);
	~Cell();
    
    //Getters
    int get_cell_number();
    Vec get_origin_of_cell();
    int get_number_of_atoms_in_cell();
    std::vector<Atom*> get_atoms_in_cell();
    
    //Setters
    void set_cell_number(int);
    void set_origin_of_cell(Vec);
    
    //Other functions
	void add_atom(Atom*);
	void clear_cell();
};

#endif