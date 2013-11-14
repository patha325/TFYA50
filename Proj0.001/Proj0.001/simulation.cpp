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
	std::string crystal_structure	: Name of chrystal structure (scc, fcc, hcp, bcc...)
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
	//Vec prev_acceleration = Vec(0,0,0); //Används ej
	
	k_b = 8.617342e-5f; //[eV][K]^{-1}
	initial_velocity_modulus = sqrt((3*k_b*temperature)/(mass));
	cout << "initial_velocity_modulus " << initial_velocity_modulus << endl;
    
    //Initial setup
    create_list_of_atoms();
	cout << "Total number of atoms: " << list_of_atoms.size() << endl;

	create_cell_list();

	number_of_atoms = list_of_atoms.size();	

	//Clear files that will be written to for every simulation.
	std::ofstream fs("atoms.txt", ios::trunc);	
	std::ofstream fs2("energytemp.txt", ios::trunc);

	// Write atom position to a file so that they can be plotted in matlab using plotter.m from drive.
	for(string::size_type i = 0; i < list_of_atoms.size();i++){
		// string::size_type ist för int eftersom .size() returnerar en unsigned int, blir varning annars.

		//cout << i<<endl;
		//cout << list_of_atoms[i]->get_position()<<endl;
		//	ofstream myfile;
		//myfile.open ("example.txt");
		std::ofstream fs("atoms.txt", ios::app); 
		fs << list_of_atoms[i]->get_position()<<endl;
		fs.close();
	}

	// Write out steps, time_step and dummy index to energytemp.
	fs2 << steps << " " << time_step << " " << 0  << " " << 0 <<endl;

	fs2.close();
		   	
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

	cout << "--------------------------------- Init -----" << endl;
	cout << "temperature " << temperature << endl;
	for(int i = 0; i < number_of_atoms; i++){
		Atom* atom = list_of_atoms[i];
		cout << "atom " << i << endl;
		cout << "    position " << atom->get_position() << endl;
		cout << "    velocity " << atom->get_velocity() << endl;
		cout << "    abs_value of velocity " << atom->get_velocity().length() << endl;
		cout << "    acceleration " << atom->get_acceleration() << endl;

		}
	cout << "-------------------------------------------" << endl << endl;

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
vector list_of_atoms. Uses help functions xcc_structure() and xcc_structure_x()
-----------------------------*/
void Simulation::create_list_of_atoms(){
	if(crystal_structure == "scc"){
		scc_structure();
	}
	else if(crystal_structure == "fcc"){
		fcc_structure();
	}
	else if(crystal_structure == "bcc"){
		bcc_structure();
	}
}
void Simulation::scc_structure(){
		for(int k=0;k<unit_cells_z;k++){//Create the cells in z
		for(int j=0;j<unit_cells_y;j++){//Create the cells in y
			scc_structure_x(j,k);// Create the cells in x
		}
	}
}
void Simulation::scc_structure_x(int j, int k)
{
	for(int i=0;i<unit_cells_x;i++){
			Vec origin (i*lattice_constant,j*lattice_constant,k*lattice_constant);
			Vec extra (0,0,0);
			Vec acceleration (0,0,0);


			float cutoff = 0.5f; // The cutoff given to all of the atoms SHOULD BE CHANGED!
			list_of_atoms.push_back(new Atom(origin,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus));


	}
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
			list_of_atoms.push_back(new Atom(origin,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma, epsilon, mass,time_step,initial_velocity_modulus));
			extra = Vec(0.5f*lattice_constant,0.5f*lattice_constant,0);
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon, mass,time_step,initial_velocity_modulus));
			extra = Vec(0,0.5f*lattice_constant,0.5f*lattice_constant);
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon, mass,time_step,initial_velocity_modulus));
			extra = Vec(0.5f*lattice_constant,0,0.5f*lattice_constant);


			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus));
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

		list_of_atoms.push_back(new Atom(origin,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus));	
		extra = Vec(0.5f*lattice_constant,0.5f*lattice_constant,0.5f*lattice_constant);


		extra +=origin;
		list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus));
	}


}

/*----------------------------
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
	
	Not every time step:
	Clear cells/cell_list
	Fill cells/cell_list

	Save energies and temperature to txt.
------------------------------*/
void Simulation::next_time_step(int current_time_step){
	
	cout << "----------------------------------- t=" << current_time_step << " -----" << endl;

	float E_pot = 0;
	float E_kin = 0;
	float temperature = 0;
	float tmp_E_kin = 0;
		
	for(int i = 0; i < number_of_atoms; i++){
		//cout << "Before timestep - atom " << i << endl;
		Atom* atom = list_of_atoms[i];
		vector<Atom*> neighbouring_atoms = cell_list->get_neighbours(atom);
		//Calculate potential energy
		E_pot += atom->calculate_potential(neighbouring_atoms);

		//cout << "    E_pot " << E_pot << endl;
		//Calculate kinetic energy
		tmp_E_kin = atom->calculate_kinetic_energy();
		//cout << "    tmp_E_kin " << tmp_E_kin << endl;

		E_kin += tmp_E_kin;
		
		//Calculate temperature
		temperature += atom->calculate_temperature(tmp_E_kin);

		//cout << "    temperature " << temperature << endl;


		//Calculate next position with help from velocity, previous acceleration and current acceleration
		atom->set_next_position(atom->calculate_next_position());
		Vec new_acceleration = atom->calculate_acceleration(neighbouring_atoms);
		
		//Update everything except position for the atom
		//Velocity
		if(current_time_step != 0){
			atom->set_velocity(atom->calculate_velocity());
		}

		//cout << "atom " << i << " velocity " << atom->get_velocity() << endl;

		//Previous position
		atom->set_prev_position(atom->get_position());		
		//Previous acceleration
		atom->set_prev_acceleration(atom->get_acceleration());
		//Acceleration
		atom->set_acceleration(new_acceleration);
	}
	

	temperature = temperature/number_of_atoms;
	total_energy = E_pot + E_kin;

	
	for(int i = 0; i < number_of_atoms; i++){
		Atom* atom = list_of_atoms[i];
		//Update position
		atom->set_position(atom->get_next_position());
	}
	
	if (fmod(current_time_step, 5.0) == 0){
			cell_list->clear_cells();
			cell_list->add_atoms_to_cells(list_of_atoms);
	}

	
	for(int i = 0; i < number_of_atoms; i++){
		Atom* atom = list_of_atoms[i];
		cout << "atom " << i << endl;
		cout << "    position " << atom->get_position() << endl;
		cout << "    velocity " << atom->get_velocity() << endl;
		cout << "    abs_value of velocity " << atom->get_velocity().length() << endl;
		cout << "    acceleration " << atom->get_acceleration() << endl;

	}
	
	
	
	cout << "total_energy " << total_energy << endl;
	cout << "E_pot " << E_pot << endl;
	cout << "E_kin " << E_kin << endl;
	cout << "temperature " << temperature << endl;
	//cout << "number of atoms " << number_of_atoms << endl << endl;
	cout << "-------------------------------------------" << endl << endl;

	// Write Energy & temp to a file so that they can be plotted in matlab using plotter.m from drive.
	std::ofstream fs2("energytemp.txt", ios::app);
	fs2 << total_energy << " " << E_pot << " " << E_kin << " " << temperature <<endl;
	fs2.close();


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












