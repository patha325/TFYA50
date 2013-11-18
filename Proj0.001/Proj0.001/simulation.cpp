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
                        float new_epsilon,float new_lattice_constant,string new_crystal_structure,bool new_thermostat, 
						map<string, vector<Vec>> new_last_state){
    
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
	last_state = new_last_state;
	//Vec prev_acceleration = Vec(0,0,0); //Används ej
	
	k_b = 8.617342e-5f; //[eV][K]^{-1}
	initial_velocity_modulus = sqrt((3*k_b*temperature)/(mass));
	cout << "wanted temperature " << temperature << endl;
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
map<string, vector<Vec>> Simulation::run_simulation(){
	//World is already created

	cout << "--------------------------------- Init -----" << endl;
	cout << "temperature " << temperature << endl;
	/*
	for(int i = 0; i < number_of_atoms; i++){
		Atom* atom = list_of_atoms[i];
		cout << "atom " << i << endl;
		cout << "    position " << atom->get_position() << endl;
		cout << "    velocity " << atom->get_velocity() << endl;
		cout << "    abs_value of velocity " << atom->get_velocity().length() << endl;
		cout << "    acceleration " << atom->get_acceleration() << endl;

		}
		*/
	cout << "-------------------------------------------" << endl << endl;

	//If this is a back to back simulation, we must update the system to the last state of previous simulation
	if (last_state.size() != 0){
		cout << "Updating atoms to last state of previous simulation" << endl;
		update_atoms_btb();
	}

	//Start next_time_step while loop
	//The checks in the while loop are to make sure that we save the last state
	int i = 0;
	while(i < steps){
		bool second_to_last_time_step = false;
		bool next_to_last_time_step = false;
		bool last_time_step = false;
		if(i == steps - 3){
			second_to_last_time_step = true;
		}
		else if(i == steps - 2){
			next_to_last_time_step = true;
		}
		else if(i == steps - 1){
			last_time_step = true;
		}
		next_time_step(i, second_to_last_time_step, next_to_last_time_step, last_time_step);
		i++;
	}

	return last_state;
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
void Simulation::next_time_step(int current_time_step, bool second_to_last_time_step, bool next_to_last_time_step, bool last_time_step){
	
	cout << "--------------------------------- t=" << current_time_step << " -----" << endl;

	float E_pot = 0;
	float E_kin = 0;
	if (!thermostat){
		float temperature = 0;
	}
	float tmp_E_kin = 0;
		
	for(int i = 0; i < number_of_atoms; i++){
		//cout << "Before timestep - atom " << i << endl;
		Atom* atom = list_of_atoms[i];
		vector<Atom*> neighbouring_atoms = cell_list->get_neighbours(atom);

		//Before doing anything else: Update atom_positions with current position for each atom for this time step
		if (fmod(current_time_step, 5.0) == 0){
			atom_positions[current_time_step].push_back(atom->get_position());
		}

		//Calculate potential energy
		E_pot += atom->calculate_potential(neighbouring_atoms);

		//cout << "    E_pot " << E_pot << endl;
		//Calculate kinetic energy
		tmp_E_kin = atom->calculate_kinetic_energy();
		//cout << "    tmp_E_kin " << tmp_E_kin << endl;

		E_kin += tmp_E_kin;

		//Calculate next position with help from velocity, previous acceleration and current acceleration
		atom->set_next_position(atom->calculate_next_position());
		Vec new_acceleration = atom->calculate_acceleration(neighbouring_atoms);
		
		//Calculate temperature if not first time step
		if(current_time_step != 0 || last_state.size() != 0){
			Vec new_velocity = atom->calculate_velocity();
			if (thermostat){
				//Update velocity so temperature is constant
				float new_velocity_modulus = new_velocity.length();
				float right_modulus = 0;
				if (new_velocity_modulus != 0){
					right_modulus = initial_velocity_modulus/new_velocity_modulus;
				}
				else{
					right_modulus = initial_velocity_modulus;
					new_velocity = atom->generate_random_vector();
				}
				atom->set_velocity(right_modulus*new_velocity);
			}
			else{
				//temperature += atom->calculate_temperature(tmp_E_kin);
				//Update velocity where total energy is constant
				atom->set_velocity(new_velocity);
			}
			temperature += atom->calculate_temperature(tmp_E_kin);
		}

		//Previous position
		atom->set_prev_position(atom->get_position());		
		//Previous acceleration
		atom->set_prev_acceleration(atom->get_acceleration());
		//Acceleration
		atom->set_acceleration(new_acceleration);

		// Check if last state should be saved to last_state
		//last_state = {"next_position":[...], "position":[...], "velocity":[...], "acceleration":[...]
		//				"prev_position":[...], "prev_acceleration":[...], "next_acceleration":[...]}
		if (second_to_last_time_step){
			last_state["prev_position"].push_back(atom->get_position());
			last_state["prev_acceleration"].push_back(atom->get_acceleration());
		}
		else if (next_to_last_time_step){
			last_state["position"].push_back(atom->get_position());
			last_state["velocity"].push_back(atom->get_velocity());
			last_state["acceleration"].push_back(atom->get_acceleration());
		}
		else if (last_time_step){
			last_state["next_position"].push_back(atom->get_position());
			last_state["next_acceleration"].push_back(atom->get_acceleration());
		}
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
	
	/*
	for(int i = 0; i < number_of_atoms; i++){
		Atom* atom = list_of_atoms[i];
		cout << "atom " << i << endl;
		cout << "    position " << atom->get_position() << endl;
		cout << "    velocity " << atom->get_velocity() << endl;
		cout << "    abs_value of velocity " << atom->get_velocity().length() << endl;
		cout << "    acceleration " << atom->get_acceleration() << endl;

	}
	*/

	// Write atom position to a file so that they can be plotted in matlab using plotter.m from drive.
	// Write to file every time step
	// Seperate the positions for different timesteps
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
FUNCTION update_atoms_btb
Paramteters: None
Returns: Nothing
- 
Changes starting state of system
to be same as last state of previous
simulation
------------------------------*/
void Simulation::update_atoms_btb(){
	for(int i = 0; i < number_of_atoms; i++){
		Atom* atom = list_of_atoms[i];

		//Previous
		atom->set_prev_position(last_state["prev_position"][i]);
		atom->set_prev_acceleration(last_state["prev_acceleration"][i]);
		//Current
		atom->set_position(last_state["position"][i]);
		atom->set_velocity(last_state["velocity"][i]);
		atom->set_acceleration(last_state["acceleration"][i]);
		//Next
		atom->set_next_position(last_state["next_position"][i]);
	}
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


std::vector<Atom*> Simulation::get_list_of_atoms(){
	
	return list_of_atoms;
}

int Simulation::get_number_of_atoms(){
	
	return number_of_atoms;
}

