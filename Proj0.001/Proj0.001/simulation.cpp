#include "simulation.h"
#include <math.h>
#include <time.h>
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
	
	
	std::ofstream fs2("energytemp.txt", ios::trunc);
	// Write out steps, time_step and dummy index to energytemp.
	fs2 << steps << " " << steps <<" "<<time_step<<" "<< 0<<" "<< 0<<" "<< 0<<" "<< 0<<" "<< 0<<" "<<0 <<endl;
	fs2.close();
	
	
	
	/*// Write atom position to a file so that they can be plotted in matlab using plotter.m from drive.
	for(unsigned int i = 0; i < list_of_atoms.size();i++){
		// unsigned int ist för int eftersom .size() returnerar en unsigned int, blir varning annars.

		cout << i<<endl;
		cout << list_of_atoms[i]->get_position()<<endl;
			ofstream myfile;
		myfile.open ("example.txt");
		
		std::ofstream fs("atoms.txt", ios::app); 
		fs << list_of_atoms[i]->get_position()<<endl;
		fs.close();
	
	}*/

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

	volume = unit_cells_x*unit_cells_y*unit_cells_z*pow(lattice_constant,3); 

	// Edit the first line of energytemp after that add other data to the file.
	configure_data(steps);

}
void Simulation::configure_data(int steps){
	std::vector<float> total_energy_vector;
	std::vector<float> pot_energy_vector;
	std::vector<float> kin_energy_vector;
	std::vector<float> temperature_vector;
	std::vector<float> pressure_vector;
	std::vector<float> MSD_vector;
	std::vector<float> Debye_temp_vector;
	std::vector<float> Diff_coeff_vector;
	std::vector<float> coh_e_vector;
	float in_t_energy;
	float in_e_pot;
	float in_e_kin;
	float in_temp;
	float in_pressure;
	float in_MSD;
	float in_Debye_temp;
	float in_Diff_coeff;
	float in_coh_e;
	istringstream iss;
	string line;
	ifstream in("energytemp.txt");
	int line_number = 0;

	// Read in all data from energytemp
	while(getline(in,line)){
		istringstream iss(line);
		if(!line_number==0){
			iss>>in_t_energy>>in_e_pot>>in_e_kin>>in_temp>>in_pressure>>in_MSD>>in_Debye_temp>>in_Diff_coeff>>in_coh_e;
			total_energy_vector.push_back(in_t_energy);
			pot_energy_vector.push_back(in_e_pot);
			kin_energy_vector.push_back(in_e_kin);
			temperature_vector.push_back(in_temp);
			pressure_vector.push_back(in_pressure);
			MSD_vector.push_back(in_MSD);
			Debye_temp_vector.push_back(in_Debye_temp);
			Diff_coeff_vector.push_back(in_Diff_coeff);
			coh_e_vector.push_back(in_coh_e);
		}
		line_number++;
	}	
	in.close();
	step_out+=steps;
	
	std::ofstream fs2("energytemp.txt", ios::trunc);
	fs2 << steps << " " << step_out << " " <<time_step<< " " << 0<< " " << 0<< " " << 0<< " " << 0<< " " << 0<<" "<< 0<<endl;
	fs2.close();
	
	for(unsigned int i=0;i<total_energy_vector.size();i++){
		std::ofstream fs2("energytemp.txt", ios::app);
		fs2 << total_energy_vector[i] << " " << pot_energy_vector[i] << " " << kin_energy_vector[i]<< " " << temperature_vector[i]<< " " <<pressure_vector[i]<< " " <<MSD_vector[i]<< " " <<Debye_temp_vector[i]<< " " <<Diff_coeff_vector[i] <<" "<<coh_e_vector[i]<<endl;
		fs2.close();
	}

}





/*-------------------------
DESTRUCTOR
Destroys all atoms and the 
cells and cell list.
-------------------------*/
Simulation::~Simulation (){
	for(int i = 0; i < number_of_atoms; i++){
		delete list_of_atoms[i];
	}
	delete cell_list;
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
		cout << "Specific heat " << C_v << endl;
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

	//Set atom number for each atom
	for (unsigned int i = 0; i < list_of_atoms.size(); i++){
		list_of_atoms[i]->set_atom_number(i);
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

			list_of_atoms.push_back(new Atom(origin,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));
		}
}

/*
void Simulation::scc_corrector(){ // Corrects for the missing atoms if there is no periodic condition in the z-axis
	for(int i=0;i<unit_cells_x;i++){
		for(int j=0;j<unit_cells_y;j++){
	Vec origin (i*lattice_constant,j*lattice_constant,unit_cells_z*lattice_constant);
			Vec extra (0,0,0);
			Vec acceleration (0,0,0);
			list_of_atoms.push_back(new Atom(origin,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));
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
			list_of_atoms.push_back(new Atom(origin,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma, epsilon, mass,time_step,initial_velocity_modulus,pbc_z));
			extra = Vec(0.5f*lattice_constant,0.5f*lattice_constant,0);
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon, mass,time_step,initial_velocity_modulus,pbc_z));
			extra = Vec(0,0.5f*lattice_constant,0.5f*lattice_constant);
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon, mass,time_step,initial_velocity_modulus,pbc_z));
			extra = Vec(0.5f*lattice_constant,0,0.5f*lattice_constant);
			
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));
	}
}

/*
void Simulation::fcc_corrector(){// Corrects for the missing atoms if there is no periodic condition in the z-axis
	for(int i=0;i<unit_cells_x;i++){
		for(int j=0;j<unit_cells_y;j++){
	Vec origin (i*lattice_constant,j*lattice_constant,unit_cells_z*lattice_constant);
			Vec extra (0,0,0);
			Vec acceleration (0,0,0);
			list_of_atoms.push_back(new Atom(origin,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));
			extra = Vec(0.5f*lattice_constant,0.5f*lattice_constant,0);
			extra +=origin;
			list_of_atoms.push_back(new Atom(extra,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon, mass,time_step,initial_velocity_modulus,pbc_z));
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

		list_of_atoms.push_back(new Atom(origin,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));	

		extra = Vec(0.5f*lattice_constant,0.5f*lattice_constant,0.5f*lattice_constant);
		
		extra +=origin;
		list_of_atoms.push_back(new Atom(extra,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));
	}
}

/*
void Simulation::bcc_corrector(){// Corrects for the missing atoms if there is no periodic condition in the z-axis
	for(int i=0;i<unit_cells_x;i++){
		for(int j=0;j<unit_cells_y;j++){
	Vec origin (i*lattice_constant,j*lattice_constant,unit_cells_z*lattice_constant);
			Vec extra (0,0,0);
			Vec acceleration (0,0,0);
			list_of_atoms.push_back(new Atom(origin,cutoff,unit_cells_x,unit_cells_y,unit_cells_z,lattice_constant,sigma,epsilon,mass,time_step,initial_velocity_modulus,pbc_z));

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

	clock_t t6 = clock();

	//if (fmod(current_time_step, 5.0) == 0){
	cout << "--------------------------------- t=" << current_time_step << " -----" << endl;
	//}

	if(current_time_step == 1295){
		cout << "time_step = 1295" << endl;
	}

	float E_pot = 0;
	float E_kin = 0;
	float temperature = 0;
	float tmp_E_kin = 0;
	float pressure = 0;
	float MSD = 0;
	float Debye_temp = 0;
	float tmp_diff_coeff = 0;
	float coh_e =0;
		
	//Update atoms' positions to next position
	Atom* atom;
	Vec new_acceleration;

	if (fmod(current_time_step, 5.0) == 0){
		cell_list->clear_cells();
	}
	for(int i = 0; i < number_of_atoms; i++){
		atom = list_of_atoms[i];
		//prev_position = position etc.
		atom->update_atom();
		//Calculate and set position
		atom->calculate_and_set_position();
		atom->clear_tmp_force();

		//Update cell_list every fifth time step
		if (fmod(current_time_step, 5.0) == 0){
			cell_list->add_atom_to_cells(atom);

		
		}
		
	}

	vector<Atom*> neighbouring_atoms;
	Atom* neighbouring_atom;

	float at = 0;
	float bt = 0;
	float ct = 0;
	float dt = 0;
	clock_t t1 = clock();

	for(int i = 0; i < number_of_atoms; i++){
		atom = list_of_atoms[i];

		clock_t t2 = clock();

		//cout << "next_time_step cell_list->get_neighbours(atom).size(): " << cell_list->get_neighbours(atom).size() << endl;
		atom->update_neighbour_list(cell_list->get_neighbours(atom));
		vector<Atom*> neighbouring_atoms = atom->get_atom_neighbours();

		clock_t t3 = clock();

		//neighbouring_atoms = cell_list->get_neighbours(atom); 

		//neighbouring_atoms = atom->reduce_neighbours_list(neighbouring_atoms);

		//Before doing anything else: Update atom_positions with current position for each atom for this time step
		//This is for animation
		/*
		if (fmod(current_time_step, 5.0) == 0){
			atom_positions[current_time_step].push_back(atom->get_position());
		}
		*/
		float E_pot_tmp = 0;
		//Calculate things which need to loop over neighbouring atoms
//		cout << "Total neighbours: " << neighbouring_atoms.size() << endl;
		int count = 0;
		for(unsigned int j = 0; j < neighbouring_atoms.size(); j++){
			neighbouring_atom = neighbouring_atoms[j];
			//Cacluate distance vector to this neighbouring atom
			Vec distance = atom->distance_vector(neighbouring_atom);
			float distance_length = distance.length();
			if(atom->get_atom_number() != neighbouring_atom->get_atom_number() && distance_length <= cutoff){
				count ++;
				//Calculate potential energy
				//times two because of how we loop over neighbouring atoms
				E_pot += 2*atom->calculate_potential(distance_length);
				
				//Calculate force
				Vec tmp_force = atom->calculate_force(distance, distance_length);
				//if(atom->get_atom_number() == 3 && neighbouring_atom->get_atom_number() == 62) {
				//	cout << "distance_length " << distance_length << endl;
				//}
				/*if(E_pot_tmp > 0){
					cout << "tmp_force stor, = " << tmp_force << endl;
					cout << "atom " << atom->get_atom_number() << "och atom " << neighbouring_atom->get_atom_number() << endl;
					cout << "distance_lenght = " << distance_length << endl;
				}*/
				atom->add_tmp_force(tmp_force);
				neighbouring_atom->add_tmp_force(-1*tmp_force);
				//Calculate pressure
				pressure += 2*atom->calculate_pressure(tmp_force, distance_length);
			}
		}
		//cout << "Neighbours within cutoff: " << count << endl;
		//system("pause");
		//Destruct neighbouring_atom?

		//Calculate and set acceleration = f(r)
		atom->calculate_and_set_acceleration(atom->get_tmp_force());
		//Calculate and set velocity = f(prev_vel, acc, prev_acc)
		calculate_and_set_velocity(atom);

		//Kinetic energy
		tmp_E_kin = atom->calculate_kinetic_energy();
		E_kin += tmp_E_kin;
	
		
		clock_t t8 = clock();
		clock_t t9 = clock();
		clock_t t10 = clock();

		//Calculate and set correct velocity only if not first time step
		if(current_time_step == 1){
				atom->set_initial_velocity(atom->get_velocity());
				atom->set_initial_position(atom->get_position());
		}

		//if(current_time_step < 500 || fmod(current_time_step,10.0) == 0){ // Görs alltid vid de 500 första tidsstegen, sedan var 10e tidssteg
		
			//Calculate temperature (even if thermostat. want to check to see that temperature is held constant)
			temperature += atom->calculate_temperature(tmp_E_kin); 
			// Used to calculate Diff_coeff and MSD, 1 shold be change to another number if we decide equilibrium is reached after "number" steps...

			if(equilibrium && current_time_step > 10){
				// Diffusion coefficient
				tmp_diff_coeff += atom->get_velocity().length()*atom->get_initial_velocity().length();
				// Mean square distance
				MSD += calculate_MSD(atom);
			}
		//}

		clock_t t4 = clock();
		at += t3 - t2;
		bt += t4 - t3;
		ct += t9 - t8;
		dt += t10 - t9;
	}
	clock_t t5 = clock();


	//if(current_time_step < 500 || fmod(current_time_step,10.0) == 0){ // Görs alltid vid de 500 första tidsstegen, sedan var 10e tidssteg
		//Calculate average temperature of system
		temperature = temperature/number_of_atoms;
		//Calculate total pressure of system
		pressure = number_of_atoms*k_b*temperature/volume + 1/(6*volume)*pressure;
		//Calculate total energy of system
		total_energy = E_pot + E_kin;
	
		//Calculate cohesive energy
		if(equilibrium && current_time_step > 10){
			coh_e = E_pot/number_of_atoms; // cohesive energy is the same as potential when equilibrium is reached.
			//Calculate average MSD
			MSD = MSD/number_of_atoms;
			//Calculate diffusion coeff
			tmp_diff_coeff = tmp_diff_coeff/number_of_atoms;
			Diff_coeff += time_step*(tmp_diff_coeff + prev_diff_coeff)/2;
			prev_diff_coeff = tmp_diff_coeff;
			//Calculate Debye temperature
			Debye_temp += 3*pow(hbar,2)*temperature/(atom->get_mass()*k_b*MSD);
		}
	
	// Write Energy & temp to a file so that they can be plotted in matlab using plotter.m from drive.
	std::ofstream fs2("energytemp.txt", ios::app);

		// Write atom position to a file so that they can be plotted in matlab using plotter.m from drive.
		// Write to file every time step
		// Seperate the positions for different timesteps
		
		fs2 << total_energy << " " << E_pot << " " << E_kin << " " << temperature << " " <<pressure<< " " << MSD<< " " <<Debye_temp<< " " <<Diff_coeff<<" "<<coh_e<<endl;
		fs2.close();
	//}

	/*
	clock_t t7 = clock();
	at = at/number_of_atoms;
	bt = bt/number_of_atoms;
	ct = ct/number_of_atoms;
	dt = dt/number_of_atoms;
	if(equilibrium){
		cout << "whole next_time_step " << t7 - t6 << endl;
		cout << "whole for loop " << t5 - t1 << endl;
		cout << "tmp_diff_coeff " << ct << endl;
		cout << "MSD " << bt << endl << endl;
	}
	*/

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
	atom->calculate_and_set_velocity();
	if (thermostat){
		//Update velocity so temperature is constant
		Vec new_velocity = atom->get_velocity();
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
	
	float MSD_tmp = 0;

	Vec temp_pos = atom->get_position();
	Vec equi_pos = atom->get_initial_position();
	float temp_diff = (temp_pos - equi_pos).length();
	MSD_tmp = pow(temp_diff,2);

	return MSD_tmp;
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
	for(unsigned int i = 0; i < list_of_atoms.size(); i++){
		cell_list->add_atom_to_cells(list_of_atoms[i]);
	}
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
