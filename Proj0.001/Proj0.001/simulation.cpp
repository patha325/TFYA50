#include "simulation.h"
#include <math.h>

using namespace std;


/*-----------------------
CONSTRUCTOR
Parameters:
	int unit_cells_x				: Number of unit cells in x-direction
	int unit_cells_y				: Number of unit cells in y-direction
	int unit_cells_z				: Number of unit cells in x-direction
	int time_step					: Length of one time step
	int steps						: number of time steps for the simulation
	float temperature				: Simulation temperature
	float cutoff					: Cutoff distance for potential calculations
	float mass						: Mass of an atom
	float sigma						: Material constant sigma
	float epsilon					: Material constan epsilon
	float lattice_constant			: Lattice constant
	std::string crystal_structure	: Name of chrystal structure (fcc, hcp, bcc...)
	bool thermostat					: If a thermostat is employed
-
Creates a simulation object. 
Calls constructors for all atoms and the cell list. 
------------------------*/

Simulation::Simulation (int new_unit_cells_x, int new_unit_cells_y, int new_unit_cells_z, float new_time_step,
                        int new_steps,float new_temperature,float new_cutoff,float new_mass,float new_sigma,
                        float new_epsilon,float new_lattice_constant,string new_crystal_structure,bool new_thermostat){
    
    //Save parameters
	unit_cells_x = new_unit_cells_x;
    unit_cells_y = new_unit_cells_y;
    unit_cells_z = new_unit_cells_z;
    time_step = new_time_step;
    steps = new_steps;
    temperature = new_temperature;
    cutoff = new_cutoff;
    mass = new_mass;
    sigma = new_sigma;
    epsilon = new_epsilon;
    lattice_constant = new_lattice_constant;
    crystal_structure = new_crystal_structure;
    thermostat = new_thermostat;
	Vec prev_acceleration = Vec(0,0,0);
    
    //Initial setup
    create_list_of_atoms();
	cout << "Total number of atoms: " << list_of_atoms.size() << endl;

	create_cell_list();

	number_of_atoms = list_of_atoms.size();

	// Atoms for testing
	Atom a(Vec(0,0,0),prev_acceleration,0,new_unit_cells_x,new_unit_cells_y,new_unit_cells_z,sigma,epsilon, mass,time_step);
	Atom b(Vec(1.8,0,0),prev_acceleration,0,new_unit_cells_x,new_unit_cells_y,new_unit_cells_z,sigma,epsilon, mass,time_step);
	Atom c(Vec(0,0,0),prev_acceleration,0,new_unit_cells_x,new_unit_cells_y,new_unit_cells_z,sigma,epsilon, mass,time_step);
	Atom d(Vec(0.1,0.1,0.1),prev_acceleration,0,new_unit_cells_x,new_unit_cells_y,new_unit_cells_z,sigma,epsilon, mass,time_step);

	vector<Atom*> atomer(1);
	atomer[0] = &b;
	//atomer[1] = &c;
	//atomer[2] = &d;

	//a.distance_vector(&b);
	a.calculate_force(atomer);
	//a.calculate_potential(&b);	

	/*for(int i=0;i<list_of_atoms.size();i++){
		cout << i<<endl;
		cout << list_of_atoms[i]->get_position()<<endl;
	}*/

	std::ofstream fs("example.txt", ios::trunc);
	for(int i=0;i<list_of_atoms.size();i++){
		//cout << i<<endl;
		//cout << list_of_atoms[i]->get_position()<<endl;
		//	ofstream myfile;
		//myfile.open ("example.txt");
		std::ofstream fs("example.txt", ios::app);
		fs << list_of_atoms[i]->get_position()<<endl;
	
		fs.close();
	}

    /*
	Vec test1 (1,2,3);
	%Atom p (test1,1);
	Vec test2 (3,-1,4);
	Atom m (test2,1);

	test1 += test2;

	cout << test1 << endl;
	*/
    //Cell* myCell = new Cell(9,Vec(1,4,7));
    //myCell->add_atom(&m);
    //vector<Atom*> atomsVector;
    //atomsVector.insert(atomsVector.begin(), myCell->get_atoms_in_cell().begin(), myCell->get_atoms_in_cell().end());

    //cout << atomsVector[0]->get_position()<< endl;
	
	// Todo: Save all the input!	
}

/*-------------------------
DESTRUCTOR
Destroys all atoms and the 
cells and cell list.
-------------------------*/
Simulation::~Simulation (){

}

/*----------------------------
FUNCTION: run_simulation
Paramteters: None
Returns: Nothing
-
Here is the main loop for the 
simulation. Handles time steps, 
and everything that happes during
the simulation.
----------------------------*/
void Simulation::run_simulation(){
	//World is already created
	int i = 0;
	while(i < steps){
		next_time_step(i);
		i++;
	}

	// Todo: How often shall we update cell list?
}

/*-----------------------------
FUNCTION: create_list_of_atoms
Parameters: None
Returns: Nothing
-
Create all atoms and add them to the
vector list_of_atoms.
-----------------------------*/
void Simulation::create_list_of_atoms(){

	if(crystal_structure == "fcc"){
		fcc_structure();
	}
	if(crystal_structure == "bcc"){
		bcc_structure();
	}
		/*
	// Calculate number of atoms
		Vec extra (0,0,0);
	if(crystal_structure == "fcc"){
	for(int i=1;i<=unit_cells_x;i++){
		Vec origin (0,0,0);
		list_of_atoms.push_back(new Atom(origin+extra,0.5));
		extra= Vec (lattice_constant,0,0);

	}
	
	}
     */
	//Todo: Lattice constant & crystal_structure & unit_cells_i
	//Todo: Create all atoms in this class? Convert from fcc to atom positions.
}

void Simulation::fcc_structure(){
	for(int k=0;k<unit_cells_z;k++){//Create the cells in z
		for(int j=0;j<unit_cells_y;j++){//Create the cells in y
			fcc_structure_x(j,k);// Create the cells in x
		}
	}
}
void Simulation::fcc_structure_x(int j, int k)
{
	for(int i=0;i<unit_cells_x;i++){
			Vec origin (i*lattice_constant,j*lattice_constant,k*lattice_constant);
			Vec extra (0,0,0);
			Vec acceleration (0,0,0);
			float cutoff = 0.5; // The cutoff given to all of the atoms SHOULD BE CHANGED!
			list_of_atoms.push_back(new Atom(origin,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,sigma, epsilon, mass,time_step));	
			extra = Vec(0.5*lattice_constant,0.5*lattice_constant,0);
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,sigma,epsilon, mass,time_step));
			extra = Vec(0,0.5*lattice_constant,0.5*lattice_constant);
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,sigma,epsilon, mass,time_step));
			extra = Vec(0.5*lattice_constant,0,0.5*lattice_constant);
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,sigma,epsilon,mass,time_step));
	}
}
void Simulation::bcc_structure(){
	for(int k=0;k<unit_cells_z;k++){//Create the cells in z
		for(int j=0;j<unit_cells_y;j++){//Create the cells in y
			bcc_structure_x(j,k);// Create the cells in x
		}
	}
}
void Simulation::bcc_structure_x(int j, int k)
{
	for(int i=0;i<unit_cells_x;i++){
		Vec origin (i*lattice_constant,j*lattice_constant,k*lattice_constant);
		Vec extra (0,0,0);
		Vec acceleration (0,0,0);
		float cutoff = 0.5; // The cutoff given to all of the atoms SHOULD BE CHANGED!
		list_of_atoms.push_back(new Atom(origin,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,sigma,epsilon,mass,time_step));	
		extra = Vec(0.5*lattice_constant,0.5*lattice_constant,0.5*lattice_constant);
		extra +=origin;
		list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,sigma,epsilon,mass,time_step));
	}


}

/*------------------------------
FUNCTION next_time_step
Paramteters: int
Returns: Nothing
-
Alter everything in the simulation to get to the next time step.
	Calculate potential and kinetic energy
	Calculate velocity
	Calculate acceleration
	Calculate next position

	Update atom's position, prev position, velocity, acceleration, prev acceleration
	
	{ Not every time step 
	Clear cells/cell_list
	Fill cells/cell_list
	}

	Save ... to txt. //Doesn't to this yet.
------------------------------*/
void Simulation::next_time_step(int current_time_step){
	
	float E_pot = 0;
	float E_kin = 0;
	float temperature = 0;
	
	for(int i = 0; i < number_of_atoms; i++){
		Atom* atom = list_of_atoms[i];

		vector<Atom*> neighbouring_atoms = cell_list->get_neighbours(atom);

		//Calculate potential energy
		E_pot += atom->calculate_potential(neighbouring_atoms);
		//Kinetic energy from velocity
		E_kin += atom->calculate_kinetic_energy();
		//Temperature
		temperature += atom->calculate_temperature(E_kin);
		//Calculate next position with help from velocity, previous acceleration and current acceleration
		atom->set_next_position(atom->calculate_next_position());
		Vec new_acceleration = atom->calculate_acceleration(neighbouring_atoms);
		
		//Update everything except position for the atom
		//Velocity
		atom->set_velocity(atom->calculate_velocity());
		//Previous position
		atom->set_prev_position(atom->get_position());
		
		//Previous acceleration
		atom->set_prev_acceleration(atom->get_acceleration());
		//Acceleration
		atom->set_acceleration(new_acceleration);
	}
	
	for(int i = 0; i < number_of_atoms; i++){
		Atom* atom = list_of_atoms[i];
		//Update position
		atom->set_position(atom->get_next_position());
		cout << "atom " << i << " previous position " << atom->get_prev_position() << endl;
		cout << "atom " << i << " position " << atom->get_position() << endl;
		cout << "atom " << i << " next position " << atom->get_next_position() << endl;
	}
	
	if (fmod(current_time_step, 5.0) == 0){
			cell_list->clear_cells();
			cell_list->add_atoms_to_cells(list_of_atoms);
	}
	
	
	temperature = temperature/number_of_atoms;
	cout << "E_pot " << E_pot << endl;
	cout << "E_kin " << E_kin << endl;
	cout << "temperature " << temperature << endl;
	cout << "number of atoms " << number_of_atoms << endl << endl;
	
	return;
}

/*------------------------------
FUNCTION regulate_thermostat
Paramteters: None
Returns: Nothing
- 
Regulates systems to maintain the 
temperature of the heat sink.
------------------------------*/
void Simulation::regulate_thermostat(){}

/*----------------------------
FUNCTION save
Parameters: None
Return: Nothing
-
Save ??? to a .txt file with some structure. 
TODO: How save, format?
----------------------------*/
void Simulation::save(){}
	
/*------------------------------
FUNCTION create_cell_list
Parameters: None
Returns: Nothing
-
Creates the cell list by calling its 
constructor. Adds all atoms to the 
cell list.
------------------------------*/
void Simulation::create_cell_list(){
	cell_list = new Cell_list(cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant);
	cell_list->add_atoms_to_cells(list_of_atoms);
}












