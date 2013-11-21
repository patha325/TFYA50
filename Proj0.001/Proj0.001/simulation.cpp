#include "simulation.h"
#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;

int step_out=0;
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
						bool new_equilibrium, bool new_pbc_z){
    
    //Save parameters
	unit_cells_x = new_unit_cells_x;
    unit_cells_y = new_unit_cells_y;
    unit_cells_z = new_unit_cells_z;
    time_step = new_time_step;
    steps = new_steps;
    temperature = new_temperature;
	pressure = 0;
    cutoff = new_cutoff;
    mass = new_mass;
    sigma = new_sigma;
    epsilon = new_epsilon;
    lattice_constant = new_lattice_constant;
    crystal_structure = new_crystal_structure;
    thermostat = new_thermostat;
	equilibrium = new_equilibrium;
	pbc_z = new_pbc_z;
	volume = unit_cells_x*unit_cells_y*unit_cells_z*pow(lattice_constant,3);
	prev_diff_coeff = 0;
	Diff_coeff = 0;
	
	//Vec prev_acceleration = Vec(0,0,0); //Används ej
	
	k_b = 8.617342e-5f; //[eV][K]^{-1}
	hbar = 0.65821189f; // [eV][fs]
	initial_velocity_modulus = sqrt((3*k_b*temperature)/(mass));
	cout << "wanted temperature " << temperature << endl;
	cout << "initial_velocity_modulus " << initial_velocity_modulus << endl;
    
    //Initial setup
    create_list_of_atoms();
	cout << "Total number of atoms: " << list_of_atoms.size() << endl;
	if (pbc_z){
		cout << "PBC is on in Z-direction." << endl;
	}
	else{
		cout << "PBC is off in Z-direction." << endl;
	} 
	cout << "Cutoff distance: " << cutoff << endl;

	create_cell_list();

	number_of_atoms = list_of_atoms.size();	

	//Clear files that will be written to for every simulation.
	//std::ofstream fs("atoms.txt", ios::trunc);
	
	//if(last_state.empty()){
	
	std::ofstream fs2("energytemp.txt", ios::trunc);
	// Write out steps, time_step and dummy index to energytemp.
	fs2 << steps << " " << steps <<" "<<time_step<<" "<< 0 <<endl;
	fs2.close();
	
	//}
	
	
	// Write atom position to a file so that they can be plotted in matlab using plotter.m from drive.
	for(string::size_type i = 0; i < list_of_atoms.size();i++){
		// string::size_type ist för int eftersom .size() returnerar en unsigned int, blir varning annars.

		//cout << i<<endl;
		//cout << list_of_atoms[i]->get_position()<<endl;
		//	ofstream myfile;
		//myfile.open ("example.txt");
		/*
		std::ofstream fs("atoms.txt", ios::app); 
		fs << list_of_atoms[i]->get_position()<<endl;
		fs.close();
	*/
	}

	step_out +=steps;
	return;  	
	// Todo: Save all the input!	
}

Simulation::Simulation(Simulation* old_simulation, int new_steps, bool new_equilibrium){

	list_of_atoms = old_simulation->get_list_of_atoms();
	number_of_atoms = old_simulation->get_number_of_atoms();
	time_step = old_simulation->get_time_step();
	steps = new_steps;
	temperature = old_simulation->get_temperature();
	cutoff = old_simulation->get_cutoff();
	thermostat = old_simulation->get_thermostat();
	pbc_z = old_simulation->get_pbc_z();
	initial_velocity_modulus = old_simulation->get_initial_velocity_modulus();
	unit_cells_x = old_simulation->get_unit_cells_x();
	unit_cells_y = old_simulation->get_unit_cells_y();
	unit_cells_z = old_simulation->get_unit_cells_z();
	total_energy = old_simulation->get_total_energy();
	cell_list = old_simulation->get_cell_list();
	equilibrium = new_equilibrium;
	prev_diff_coeff = 0;
	Diff_coeff = 0;

	//Boltzmann constant
	k_b = 8.617342e-5f;

	//Material
	mass = old_simulation->get_mass();
	sigma = old_simulation->get_sigma();
	epsilon = old_simulation->get_epsilon();
	lattice_constant = old_simulation->get_lattice_constant();
	crystal_structure = old_simulation->get_crystal_structure();


	// Edit the first line of energytemp after that add other data to the file.
	configure_data(steps);

}
void Simulation::configure_data(int steps){
	std::vector<float> total_energy_vector;
	std::vector<float> pot_energy_vector;
	std::vector<float> kin_energy_vector;
	std::vector<float> temperature_vector;
	float in_t_energy;
	float in_e_pot;
	float in_e_kin;
	float in_temp;
	istringstream iss;
	string line;
	ifstream in("energytemp.txt");
	int line_number = 0;

	// Read in all data from energytemp
	while(getline(in,line)){
		istringstream iss(line);
		if(!line_number==0){
			iss>>in_t_energy>>in_e_pot>>in_e_kin>>in_temp;
			total_energy_vector.push_back(in_t_energy);
			pot_energy_vector.push_back(in_e_pot);
			kin_energy_vector.push_back(in_e_kin);
			temperature_vector.push_back(in_temp);
		}
		line_number++;
	}	
	in.close();
	step_out+=steps;

	std::ofstream fs2("energytemp.txt", ios::trunc);
	fs2 << steps << " " << step_out << " " <<time_step<< " " << 0 <<endl;
	fs2.close();
	
	for(unsigned int i=0;i<total_energy_vector.size();i++){
		std::ofstream fs2("energytemp.txt", ios::app);
		fs2 << total_energy_vector[i] << " " << pot_energy_vector[i] << " " << kin_energy_vector[i]<< " " << temperature_vector[i] <<endl;
		fs2.close();
	}

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

	//If this is a back to back simulation, we must update the system to the last state of previous simulation
	/*
	if (last_state.size() != 0){
		update_atoms_btb();
	}
	*/

	//Start next_time_step while loop

	int i = 0;
	while(i < steps){
		next_time_step(i);
		i++;
	}
	

	// Calculate specific heat coeff
	float C_v;
	if(equilibrium) {
		C_v = calculate_specific_heat();
	}

	//Check if system is in equilibrium
	//check_equilibrium();
	return;
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
		/*
		if(!pbc_z)
			scc_corrector();*/
	}
	else if(crystal_structure == "fcc"){
		fcc_structure();
		/*
		if(!pbc_z)
			fcc_corrector();
			*/
	}
	else if(crystal_structure == "bcc"){
		bcc_structure();
		/*
		if(!pbc_z)
			bcc_corrector();*/
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

			list_of_atoms.push_back(new Atom(origin,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));
		}
}

/*
void Simulation::scc_corrector(){ // Corrects for the missing atoms if there is no periodic condition in the z-axis
	for(int i=0;i<unit_cells_x;i++){
		for(int j=0;j<unit_cells_y;j++){
	Vec origin (i*lattice_constant,j*lattice_constant,unit_cells_z*lattice_constant);
			Vec extra (0,0,0);
			Vec acceleration (0,0,0);
			list_of_atoms.push_back(new Atom(origin,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));
		}
	}
}
*/

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
			list_of_atoms.push_back(new Atom(origin,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma, epsilon, mass,time_step,initial_velocity_modulus,pbc_z));
			extra = Vec(0.5f*lattice_constant,0.5f*lattice_constant,0);
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon, mass,time_step,initial_velocity_modulus,pbc_z));
			extra = Vec(0,0.5f*lattice_constant,0.5f*lattice_constant);
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon, mass,time_step,initial_velocity_modulus,pbc_z));
			extra = Vec(0.5f*lattice_constant,0,0.5f*lattice_constant);
			
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));
	}
}

/*
void Simulation::fcc_corrector(){// Corrects for the missing atoms if there is no periodic condition in the z-axis
	for(int i=0;i<unit_cells_x;i++){
		for(int j=0;j<unit_cells_y;j++){
	Vec origin (i*lattice_constant,j*lattice_constant,unit_cells_z*lattice_constant);
			Vec extra (0,0,0);
			Vec acceleration (0,0,0);
			list_of_atoms.push_back(new Atom(origin,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));
			extra = Vec(0.5f*lattice_constant,0.5f*lattice_constant,0);
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon, mass,time_step,initial_velocity_modulus,pbc_z));
		}
	}
}
*/

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

		list_of_atoms.push_back(new Atom(origin,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));	

		extra = Vec(0.5f*lattice_constant,0.5f*lattice_constant,0.5f*lattice_constant);
		
		extra +=origin;
		list_of_atoms.push_back(new Atom(extra,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));
	}
}

/*
void Simulation::bcc_corrector(){// Corrects for the missing atoms if there is no periodic condition in the z-axis
	for(int i=0;i<unit_cells_x;i++){
		for(int j=0;j<unit_cells_y;j++){
	Vec origin (i*lattice_constant,j*lattice_constant,unit_cells_z*lattice_constant);
			Vec extra (0,0,0);
			Vec acceleration (0,0,0);
			list_of_atoms.push_back(new Atom(origin,acceleration,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));

		}
	}
}
*/

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

	if (fmod(current_time_step, 5.0) == 0){
	cout << "--------------------------------- t=" << current_time_step << " -----" << endl;
	}

	float E_pot = 0;
	float E_kin = 0;
	float temperature = 0;
	float tmp_E_kin = 0;
	float pressure = 0;
	float MSD = 0;
	float Debye_temp = 0;
	float tmp_diff_coeff = 0;
		
	//Update atoms' positions to next position
	Atom* atom;
	Vec new_acceleration;

	for(int i = 0; i < number_of_atoms; i++){
		atom = list_of_atoms[i];
		//Update position
		atom->set_position(atom->get_next_position());
	}

	vector<Atom*> neighbouring_atoms;
	Atom* neighbouring_atom;

	for(int i = 0; i < number_of_atoms; i++){
		atom = list_of_atoms[i];

		neighbouring_atoms = cell_list->get_neighbours(atom); 

		//neighbouring_atoms = atom->reduce_neighbours_list(neighbouring_atoms);

		//Before doing anything else: Update atom_positions with current position for each atom for this time step
		//This is for animation
		/*
		if (fmod(current_time_step, 5.0) == 0){
			atom_positions[current_time_step].push_back(atom->get_position());
		}
		*/

		Vec force = Vec (0, 0, 0);
		//Calculate things which need to loop over neighbouring atoms
//		cout << "Total neighbours: " << neighbouring_atoms.size() << endl;
		int count = 0;
		for(unsigned int j = 0; j < neighbouring_atoms.size(); j++){
			neighbouring_atom = neighbouring_atoms[j];
			//Cacluate distance vector to this neighbouring atom
			Vec distance = atom->distance_vector(neighbouring_atom);
			if(distance.length() != 0 && distance.length() <= cutoff){
				count ++;
				//Calculate potential energy
				E_pot += atom->calculate_potential(neighbouring_atom, distance);
				//Calculate force
				Vec tmp_force = atom->calculate_force(neighbouring_atom, distance);
				force += tmp_force;
				//Calculate pressure
				pressure += atom->calculate_pressure(neighbouring_atom, tmp_force, distance);
			}
		}
//		cout << "Neighbours within cutoff: " << count << endl;
//		system("pause");
		//Destruct neighbouring_atom?

		//Calculate new_acceleration = f(r)
		new_acceleration = atom->calculate_acceleration(force);
		//Kinetic energy
		tmp_E_kin = atom->calculate_kinetic_energy();
		E_kin += tmp_E_kin;
		//Next position
		atom->set_next_position(atom->calculate_next_position());
		
		//Calculate and set correct velocity only if not first time step
		if(current_time_step != 0){
			calculate_and_set_velocity(atom);

			if(current_time_step == 1) {
				atom->set_initial_velocity(atom->get_velocity());
			}

			if(equilibrium) {
				// Diffusion coefficient, later??
				tmp_diff_coeff += atom->get_velocity().length()*atom->get_initial_velocity().length();
				// Mean square distance
				MSD += calculate_MSD(atom);
			}
		}

		//Calculate temperature (even if thermostat. want to check to see that temperature is held constant)
		temperature += atom->calculate_temperature(tmp_E_kin);

		//Previous position
		atom->set_prev_position(atom->get_position());
		//Previous acceleration
		atom->set_prev_acceleration(atom->get_acceleration());
		//Acceleration
		atom->set_acceleration(new_acceleration);
	}

	//Calculate average temperature of system
	temperature = temperature/number_of_atoms;
	//Calculate average MSD
	MSD = MSD/number_of_atoms;
	//Calculate diffusion coeff
	tmp_diff_coeff = tmp_diff_coeff/number_of_atoms;
	Diff_coeff += time_step*(tmp_diff_coeff + prev_diff_coeff)/2;
	prev_diff_coeff = tmp_diff_coeff;
	//Calculate Debye temperature
	Debye_temp += 3*pow(hbar,2)*temperature/(atom->get_mass()*k_b*MSD);
	//Calculate total pressure of system
	pressure = number_of_atoms*k_b*temperature/volume + 1/(6*volume)*pressure;
	//Calculate total energy of system
	total_energy = E_pot + E_kin;
	
	//Calculate cohesive energy
	if(equilibrium){
		float Diff_coef = 0;
		float coh_e = E_pot; // cohesive energy is the same as potential when equilibrium is reached.
	}
	
	//Update cell_list every fifth time step
	if (fmod(current_time_step, 5.0) == 0){
			cell_list->clear_cells();
			cell_list->add_atoms_to_cells(list_of_atoms);
	}

	// Write atom position to a file so that they can be plotted in matlab using plotter.m from drive.
	// Write to file every time step
	// Seperate the positions for different timesteps
	
	// Write Energy & temp to a file so that they can be plotted in matlab using plotter.m from drive.
	//TODO: MSD and Debye_temp should also be written to file

	std::ofstream fs2("energytemp.txt", ios::app);
	
	fs2 << total_energy << " " << E_pot << " " << E_kin << " " << temperature <<endl;
	fs2.close();

	return;
}

/*------------------------------
FUNCTION calculate_and_set_velocity
Paramteters: None
Returns: None
- 
Calculates and sets correct velocity
Checks if system has thermostat
------------------------------*/
void Simulation::calculate_and_set_velocity(Atom* atom){
	//Calculate velocity as if total energy is to be constant
	Vec new_velocity = atom->calculate_velocity();
	if (thermostat){
		//Update velocity so temperature is constant
		float new_velocity_modulus = new_velocity.length();
		float right_modulus = 0;
		//This is just a simple check in order to not get infinities
		if (new_velocity_modulus != 0){
			//initial_velocity_modulus is calculated from initial temperature
			right_modulus = initial_velocity_modulus/new_velocity_modulus;
		}
		else{
			right_modulus = initial_velocity_modulus;
			//Velocity of any atom will never be zero. Is this a problem?
			new_velocity = atom->generate_random_vector();
		}
		atom->set_velocity(right_modulus*new_velocity);
	}
	else{ //not thermostat
		//Update velocity where total energy is constant
		atom->set_velocity(new_velocity);
	}
}

/*------------------------------
FUNCTION calculate_specific_heat
Paramteters: None
Returns: float specific heat
- 
Calculates soecific heat
------------------------------*/
float Simulation::calculate_specific_heat(){
	std::istringstream iss;
	string line;
	ifstream in("energytemp.txt");
	float temp=0;
	float temp_2=0;
	float tmp_temp;
	float dump;
	int number_of_time_steps = 0;
	while(getline(in,line)){
		istringstream iss(line);
		iss>>dump>>dump>>dump>>tmp_temp;
		temp+=tmp_temp;
		temp_2+=pow(tmp_temp,2);
		number_of_time_steps++;
	}
	in.close();

	float temp_av = temp/number_of_time_steps;
	float temp_2_av = temp_2/number_of_time_steps;

	float C_v;
	C_v = (3*number_of_atoms*k_b)/2*(1/(1-(temp_2_av - pow(temp_av,2))/pow(temp_av,2)*2*number_of_atoms/3));
	return C_v;
}

/*------------------------------
FUNCTION calculate_MSD
Paramteters: Atom*
Returns: float mean square displacement
- 
Calculates mean square displacement for an atom
------------------------------*/
float Simulation::calculate_MSD(Atom* atom){
	
	float MSD = 0;
	for(int i = 1; i <= number_of_atoms; i++) {
		Vec temp_pos = atom->get_position();
		Vec equi_pos (0,0,0); // Här ska positionen då jämvikt uppnåddes hämtas
		float temp_diff = (temp_pos - equi_pos).length();
		MSD += 1/number_of_atoms*pow(temp_diff,2); 
	}
	return MSD;
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
	/*
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
	*/
}

/*------------------------------
FUNCTION update_last_state()
Paramteters: None
Returns: Nothing
- 
Saves data to last state for btb simulation
------------------------------*/
void Simulation::update_last_state(Atom* atom, bool second_to_last_time_step, bool next_to_last_time_step, bool last_time_step){
	/*
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
	*/
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
	cell_list = new Cell_list(cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,pbc_z);
	cell_list->add_atoms_to_cells(list_of_atoms);
}

/*-----
GETTERS
-----*/
vector<Atom*> Simulation::get_list_of_atoms(){return list_of_atoms;}
int Simulation::get_number_of_atoms(){return number_of_atoms;}
float Simulation::get_time_step(){return time_step;}
int Simulation::get_steps(){return steps;}
float Simulation::get_temperature(){return temperature;}
float Simulation::get_cutoff(){return cutoff;}
bool Simulation::get_thermostat(){return thermostat;}
bool Simulation::get_pbc_z(){return pbc_z;}
float Simulation::get_initial_velocity_modulus(){return initial_velocity_modulus;}
int Simulation::get_unit_cells_x(){return unit_cells_x;}
int Simulation::get_unit_cells_y(){return unit_cells_y;}
int Simulation::get_unit_cells_z(){return unit_cells_z;}
float Simulation::get_total_energy(){return total_energy;}
float Simulation::get_mass(){return mass;}
float Simulation::get_sigma(){return sigma;}
float Simulation::get_epsilon(){return epsilon;}
float Simulation::get_lattice_constant(){return lattice_constant;}
string Simulation::get_crystal_structure(){return crystal_structure;}
Cell_list* Simulation::get_cell_list(){return cell_list;}
