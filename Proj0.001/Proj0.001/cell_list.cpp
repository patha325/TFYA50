#include <iostream>
#include <time.h>

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
Create cells! through cell.cpp.
------------------------------ */
Cell_list::Cell_list(float new_cutoff, int unit_cells_x, int unit_cells_y, int unit_cells_z, float new_lattice_constant, bool new_pbc_z){

    cutoff = new_cutoff; 
    lattice_constant = new_lattice_constant;
	pbc_z = new_pbc_z;
	number_of_cells = 0;
    
	/*
    cout << "Cutoff: " << cutoff << endl;
    cout << "Lattice constant: " << lattice_constant << endl;
	if(pbc_z) cout << "Periodic boundary condition IS used." << endl << endl;
	else cout << "Periodic boundary condition is NOT used." << endl << endl;
      
    cout << "Unit cells in X: " << unit_cells_x << endl;
    cout << "Unit cells in Y: " << unit_cells_y << endl;
    cout << "Unit cells in Z: " << unit_cells_z << endl << endl;
	*/

    
    bulk_length_x = unit_cells_x*lattice_constant;
    bulk_length_y = unit_cells_y*lattice_constant;
    bulk_length_z = unit_cells_z*lattice_constant;
    

    cout << "Bulk length X: " << bulk_length_x << endl;
    cout << "Bulk length Y: " << bulk_length_y << endl;
    cout << "Bulk length Z: " << bulk_length_z << endl << endl;

    
	//Could be an issue if cell_length:s are a multiple of lattice_constant
    cell_length_x = bulk_length_x/floor(bulk_length_x/cutoff);
    cell_length_y = bulk_length_y/floor(bulk_length_y/cutoff);
    cell_length_z = bulk_length_z/floor(bulk_length_z/cutoff);
    
	/*
    cout << "Cell length X: " << cell_length_x << endl;
    cout << "Cell length Y: " << cell_length_y << endl;
    cout << "Cell length Z: " << cell_length_z << endl << endl;

	system("pause");
	*/
    create_cells();
    
}


/* ------------------------------
DESTRUCTOR
PARAMETERS: None
-
Destroys cells in list_of_cells
------------------------------ */
Cell_list::~Cell_list(){

    for (unsigned int i = 0; i<list_of_cells.size(); i++) {
        delete list_of_cells[i];
    }
}

/* ---
PUBLIC
--- */

/* ------------------------------
FUNCTION: Cell_list::add_atom_to_cells()
PARAMETERS: vector<Atoms*>
RETURN: void
-
Adds all Atoms to the cells
------------------------------ */
void Cell_list::add_atom_to_cells(Atom* current_atom){


/*
    for (int i = 0 ; i<atoms_list.size(); i++) {
        cout << "Atom " << i << ": " << atoms_list[i]->get_position() << endl;
    }


    for (int i = 0; i<list_of_cells.size(); i++) {
        cout << "Cell " << i << ": "<< list_of_cells[i]->get_origin_of_cell() << endl;
    }
*/

    /*
	for (unsigned int i = 0; i< atoms_list.size(); i++) {
        Atom* current_atom = atoms_list[i];
	bool skogen = false;
	if (current_atom->get_position().getX() > bulk_length_x ||
		current_atom->get_position().getY() > bulk_length_y ||
		current_atom->get_position().getZ() > bulk_length_z) {
			cout << "ÅT SKOGEN" << endl;
			skogen = true;
	}*/

	/*Vec tmp_position = current_atom->get_position();
	//Check if atom outside cells
	//x
	while(tmp_position.getX() < 0 || tmp_position.getX() > bulk_length_x){
		int sign = my_sign(tmp_position.getX());
		tmp_position.setX(tmp_position.getX()-sign*bulk_length_x);
	}
	//y
	while(tmp_position.getY() < 0 || tmp_position.getY() > bulk_length_y){
		int sign = my_sign(tmp_position.getY());
		tmp_position.setY(tmp_position.getY()-sign*bulk_length_y);
	}
	//z
	while(tmp_position.getZ() < 0 || tmp_position.getZ() > bulk_length_z){
		int sign = my_sign(tmp_position.getZ());
		tmp_position.setZ(tmp_position.getZ()-sign*bulk_length_z);
	}*/

	/*if(skogen){
		cout << "tmp_position " << tmp_position << endl;
	}*/

	//cout << "Before while loop in add_atom_to_cells" << endl;
	int cell_number_iterator = 0;
    bool found = false;
    while (!found) {
        if (current_atom->get_position().getX()<list_of_cells[cell_number_iterator]->get_origin_of_cell().getX()+cell_length_x && 
            current_atom->get_position().getY()<list_of_cells[cell_number_iterator]->get_origin_of_cell().getY()+cell_length_y &&
            current_atom->get_position().getZ()<list_of_cells[cell_number_iterator]->get_origin_of_cell().getZ()+cell_length_z) {

                
                
//				cout << "Atom " << current_atom->get_atom_number() << " with origin " << current_atom->get_position() << " is in Cell with origin " << list_of_cells[cell_number_iterator]->get_origin_of_cell() << ": " << cell_number_iterator << endl;
                
				list_of_cells[cell_number_iterator]->add_atom(current_atom);
				//list_of_cells[cell_number_iterator]->add_atom_number(current_atom->get_atom_number());


//				cout << "Cell with number " << cell_number_iterator << " has " << list_of_cells[cell_number_iterator]->get_number_of_atoms_in_cell() << " atoms in it." << endl;

            found = true;
        }
        else {
            cell_number_iterator++;
        }
    }
	//cout << "After while loop in add_atom_to_cells" << endl;
	current_atom->set_cell_number(cell_number_iterator);
}

/* ------------------------------
FUNCTION: sign
PARAMETERS: float
RETURN: int
-
Returns +1 or -1 depending on the
sign of the incoming float
------------------------------
int Cell_list::my_sign(float number){
	int sign;
	if(number < 0){
		sign = -1;
	}
	else{
		sign = +1;
	}
	return sign;
}*/

/* ------------------------------
FUNCTION: Cell_list::get_neighbours()
PARAMETERS: Atom* atom
RETURN: vector<Atom*>
-
Returns a vector with all Atoms which
are neighbours to the parameter Atom.
------------------------------ */
vector<Atom*> Cell_list::get_neighbours(Atom* atom){
	//clock_t t1 = clock();
	//cout << "Atom number " << atom->get_atom_number() << endl;

    int cell_number = atom->get_cell_number();
	int atom_number = atom->get_atom_number();
    vector<Cell*> neighbouring_cells = number_to_cell_vector_map[cell_number];
	vector<Atom*> neighbouring_atoms;
		

	//clock_t t3 = clock();
    for (unsigned int i = 0; i < neighbouring_cells.size(); i++) {
		vector<Atom*> atoms_to_add = neighbouring_cells[i]->get_atoms_in_cell();
		for (unsigned int j = 0; j < atoms_to_add.size(); j++){
			//float distance = atom->distance_vector(atoms_to_add[j]).length();
			if (atom_number < atoms_to_add[j]->get_atom_number()){
				neighbouring_atoms.push_back(atoms_to_add[j]);
			}
		}
    }
	//clock_t t2 = clock();
	//cout << "before for loop " << t3 - t1 << endl;
	//cout << "Time for getting neighbours: " << t2 - t1 << endl;
	return neighbouring_atoms;
	
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
void Cell_list::clear_cells(){

	for(unsigned int i=0; i< list_of_cells.size(); i++){
	list_of_cells[i]->clear_cell();	 
	}
}

/* ------------------------------
FUNCTION: Cell_list::create_cells()
PARAMETERS: None
RETURN: void
-
Creates Cell-objects and adds them
to the list_of_cells vector.
 
Also creates a mapping between cell
number and a vector containing all
neighbouring cells. (number_to_cell_vector_map)
------------------------------ */
void Cell_list::create_cells(){

    Vec current_origin(0,0,0);
    Vec current_matrix_coordinates(0,0,0);
    map<Vec, int> vec_to_number_map;
    
    int i = 0; //Cell number
    float max_orgin_x; //Highest x-value a cell-origin can take
    float max_orgin_y; //Highest y-value a cell-origin can take

    float max_orgin_z; //Highest z-value a cell-origin can take    


	/*
	cout << "Max origin x: " << max_orgin_x << endl;
	cout << "Max origin x: " << max_orgin_y << endl;
	cout << "Max origin x: " << max_orgin_z << endl;
	*/
	/*
	cout  << "Bulk x: " << bulk_length_x << endl;
	cout  << "Bulk y: " << bulk_length_y << endl;
	cout  << "Bulk z: " << bulk_length_z << endl;
	*/
	/*
	cout << "Cell length x: " << cell_length_x << endl;
	cout << "Cell length y: " << cell_length_y << endl;
	cout << "Cell length z: " << cell_length_z << endl;
    */

    
    /*----
     Creates and inserts cells into list. Gives each cell a number and 
     an origin.
        Creates a mapping between matrix coordinates and cell number 
     (only used when creating the cell list structure.
     ---*/
	

	cout << endl << "Creating cells... " << endl;

    while (current_origin.getZ()<bulk_length_z){
        while (current_origin.getY()<bulk_length_y) {
            while (current_origin.getX()<bulk_length_x) {
            
                list_of_cells.insert(list_of_cells.end(), new Cell(i,current_origin));
                vec_to_number_map.insert(pair<Vec,int>(current_matrix_coordinates,i));
				number_of_cells++;
                
                //cout << i << ": " << current_origin << endl;
            
                i++;
                current_origin = (current_origin + Vec(cell_length_x,0,0));
                current_matrix_coordinates = current_matrix_coordinates + Vec(1,0,0);
                max_orgin_x = current_origin.getX();
            }
            current_origin = (current_origin+Vec(-current_origin.getX(),cell_length_y,0));
            current_matrix_coordinates = current_matrix_coordinates + Vec(-current_matrix_coordinates.getX(),1,0);
            max_orgin_y = current_origin.getY();
        }
        current_origin = (current_origin+Vec(-current_origin.getX(),-current_origin.getY(),cell_length_z));
        current_matrix_coordinates = current_matrix_coordinates + Vec(-current_matrix_coordinates.getX(),-current_matrix_coordinates.getY(),1);
        max_orgin_z = current_origin.getZ();
    }


	cout << "Cells created!" << endl;
	cout << "Created in total " << i << " cells." << endl;
	cout << "Cell dimensions are: " << cell_length_x << " x " << cell_length_y << " x "  << cell_length_z << endl << endl;

    
    //Number of cells in each direction
    int cells_x = int(max_orgin_x/cell_length_x);
    int cells_y = int(max_orgin_y/cell_length_y);
    int cells_z = int(max_orgin_z/cell_length_z);
 
/*
    cout << "Numbers of cells X: " << cells_x << endl;
    cout << "Numbers of cells Y: " << cells_y << endl;
    cout << "Numbers of cells Z: " << cells_z << endl << endl;
*/
    
    /*---
     Creates a mapping from cell number to a vector of all neighbouring cells.
     ---*/
    for (float z = 0; z<cells_z; z++){
        for (float y = 0; y<cells_y; y++){
            for (float x = 0; x<cells_x; x++) {
                float minus_x = x-1;  // Dessa kan bli int om Vec inte används nedan...
                float minus_y = y-1;
                float minus_z = z-1;
                float plus_x = x+1;
                float plus_y = y+1;
                float plus_z = z+1;
                if (minus_x<0) minus_x += cells_x;
                if (minus_y<0) minus_y += cells_y;
                if (minus_z<0) minus_z += cells_z;
                if (plus_x>=cells_x) plus_x = 0;
                if (plus_y>=cells_y) plus_y = 0;
                if (plus_z>=cells_z) plus_z = 0;
                
                int cell_number = vec_to_number_map[Vec(x,y,z)];
                vector<Cell*> cell_vector;
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(minus_x   ,minus_y    ,minus_z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(x         ,minus_y    ,minus_z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(plus_x    ,minus_y    ,minus_z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(minus_x   ,y          ,minus_z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(x         ,y          ,minus_z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(plus_x    ,y          ,minus_z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(minus_x   ,plus_y     ,minus_z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(x         ,plus_y     ,minus_z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(plus_x    ,plus_y     ,minus_z)]]);
                
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(minus_x   ,minus_y    ,z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(x         ,minus_y    ,z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(plus_x    ,minus_y    ,z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(minus_x   ,y          ,z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(x         ,y          ,z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(plus_x    ,y          ,z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(minus_x   ,plus_y     ,z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(x         ,plus_y     ,z)]]);
                cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(plus_x    ,plus_y     ,z)]]);
                
				if(!(!pbc_z && plus_z==0)){
					cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(minus_x   ,minus_y    ,plus_z)]]);
					cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(x         ,minus_y    ,plus_z)]]);
					cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(plus_x    ,minus_y    ,plus_z)]]);
					cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(minus_x   ,y          ,plus_z)]]);
					cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(x         ,y          ,plus_z)]]);
					cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(plus_x    ,y          ,plus_z)]]);
					cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(minus_x   ,plus_y     ,plus_z)]]);
					cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(x         ,plus_y     ,plus_z)]]);
					cell_vector.push_back(list_of_cells[vec_to_number_map[Vec(plus_x    ,plus_y     ,plus_z)]]);
				}
                
               number_to_cell_vector_map.insert(pair<int,vector<Cell*>>(cell_number,cell_vector));
            }
        }
    }
 
/*
	for (int cell_number = 1; cell_number < cells_x*cells_y*cells_z; cell_number++){
		cout << "Cell number: " << cell_number << ": " << endl;
		vector<Cell*> cells =  number_to_cell_vector_map[cell_number];
		for (int i = 0; i<cells.size(); i++) {
			if (i%9 ==0 && i!=0) cout << endl;
			cout << cells[i]->get_cell_number() << ", ";
		}
		cout << endl << endl;
	}
*/

}

/*
void Cell_list::print_my_cell_number(int atom_number){

	for(int i = 0; i<list_of_cells.size(); i++){
		if
	}
	cout << 
}
*/

void Cell_list::print_my_cell_neighbours(int atom_number){

	vector<Cell*> lst = number_to_cell_vector_map[atom_number];
	cout << endl << "The neighbouring cells for atom #" << atom_number << " are: ";
	for (int i = 0; i < lst.size(); i++){
		cout << lst[i]->get_cell_number() << " ";
	}
	cout << endl;
}

void Cell_list::print_atoms_in_each_cell(){

	cout << endl << "ATOMS IN EACH CELL:" << endl;
	for(int i = 0; i < number_of_cells; i++){
		cout << "Atoms in cell " << i << endl;
		vector<Atom*> list_of_atoms = list_of_cells[i]->get_atoms_in_cell();

		cout << "Atoms: ";
		for(int j = 0; j < list_of_atoms.size(); j++){
			cout << list_of_atoms[j]->get_atom_number() << " ";
		}
		cout << endl;
	
	}
}

	
Cell* Cell_list::get_cell_with_number(int cell_number){

	return list_of_cells[cell_number];
}

/*
void Cell_list::add_to_map(int key, vector<Cell*> entry){
    
    number_to_cell_vector_map.insert(pair<int, vector<Cell*>>(key,entry));
}
*/

int Cell_list::get_number_of_cells(){

	return number_of_cells;
}
