#include "simulation.h"
//#include "simulation2.h"
#include <iostream>
#include "vec.h"
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include "GraphicsTestProject.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

using namespace std;
	//initialvalues for get_input_file()
	string atom;
	string structure;
	float sigma;
	float epsilon;
	float mass;
	float lattice;
	string file = "data.txt";

void get_input_file(string intresting,string file){
istringstream iss;
string line;
ifstream in(file);
while(getline(in,line)){
	istringstream iss(line);
	iss>>atom>>sigma>>epsilon>>structure>>mass>>lattice;
	if(atom==intresting){
	break;
		}
	}
in.close();
}


int main(int argc, char** argv){
    
/*
Parameters:
	int unit_cells_x				: Number of unit cells in x-direction
	int unit_cells_y				: Number of unit cells in y-direction
	int unit_cells_z				: Number of unit cells in x-direction
	float time_step					: Length of one time step
	int steps						: number of time steps for the simulation
    float temperature				: Simulation temperature
    float cutoff					: Cutoff distance for potential calculations
    float mass						: Mass of an atom
    float sigma						: Material constant sigma LJ-potential
    float epsilon					: Material constan epsilon
    float lattice_constant			: Lattice constant
    std::string crystal_structure	: Name of chrystal structure (scc, fcc, hcp, bcc...)
    bool thermostat					: If a thermostat is employed
 */

	//initializegraphics
	float a = 0.2f; // Watch out for double warning!
	
	int input_x;
	int input_y;
	int input_z;
	float input_time_step;
	int input_steps;
	float input_temperature = a;
	float input_cutoff;
	float input_mass = a;
	float input_sigma = a;
	float input_epsilon = a;
	float input_lattice_constant = a;
	string input_crystal_structure = "fcc";
	string input_material;
	bool input_thermostat = false;
	bool input_equilibrium = false;
	bool pbc_z = false;
	
	cout << "Input the number of unit cells in x,y and z direction:" <<endl;
	cin >> input_x;
	cin >> input_y;
	cin >> input_z;
	cout << "Choose (Yes = 1/No = 0) for periodic boundry condition in z"<<endl;
	cin >> pbc_z;
	cout << "Input the wanted time step size:" << endl; 
	cin >> input_time_step;
	cout << "Input the wanted number of steps:" << endl;
	cin >> input_steps;
	cout << "Input wanted material:" << endl;
	cin >> input_material;
	get_input_file(input_material,file);
	cout << "Start temperature (K):" << endl;
	cin >> input_temperature;
	cout << "Input cutoff multiples of lattice_constant:" <<endl;
	cin >> input_cutoff;
	input_cutoff = input_cutoff*lattice;
	cout << "Simulate with thermostat? (Yes = 1/No = 0)" << endl;
	cin >> input_thermostat;


	cout << endl << "------------" << endl;
	cout << "- RUNNING -" << endl;
	cout << "-----------" << endl << endl;
	
	input_sigma = sigma; //[Å]
	input_epsilon = epsilon; //[eV]
	input_crystal_structure = structure;
	input_mass = mass; //[eV/c^2]=[eV][Å]^2[fs]^-2
	input_lattice_constant=lattice;

	cout << "Sigma: "<< sigma << endl;
	cout << "Epsilon: "<<epsilon<<endl;
	cout << "Structure: "<<structure<<endl;
	cout << "Mass: "<<mass<<endl;
	cout << "Lattice: "<<lattice<<endl;
	
		//Create first simulation world
	Simulation* simulation = new Simulation(input_x,input_y,input_z,input_time_step,input_steps,input_temperature, input_cutoff, 
											input_mass, input_sigma, input_epsilon, input_lattice_constant,input_crystal_structure,
											input_thermostat, input_equilibrium, pbc_z);
	cout << "Running simulation..." << endl << endl;
	
	simulation->run_simulation();

/*
	Simulation* simulation = new Simulation(simulation2);
	cout << "Number of atoms: " << simulation->get_list_of_atoms().size() << endl;
	cout << "Energy: " << simulation->get_total_energy() << endl;
*/

	//Run back to back simulation
	//Always with most recent simulation as it is now
	//btb = back to back

	bool back_to_back = true;
	cout << "Do you wish to run a new simulation back to back? (Yes = 1/No = 0)" << endl;
	cin >> back_to_back;
	while (back_to_back){
		cout << "System in equilibrium? (Yes = 1/No = 0)" << endl;
		cin >> input_equilibrium;
		cout << "Input the wanted number of steps:" << endl;
		cin >> input_steps;
		cout << "Starting new simulation back to back with previous!" << endl;
		//Create new simulation
		Simulation* btb_simulation = new Simulation(simulation, input_steps, input_equilibrium);
		cout << "Running simulation..." << endl << endl;
		btb_simulation->run_simulation();

		//Ask again to runt new btb simulation
		cout << "Do you wish to run a new simulation back to back? (Yes = 1/No = 0)" << endl;
		cin >> back_to_back;
	}

	system("pause");

//	plotter(argc, argv,simulation2->get_list_of_atoms(),simulation2->get_number_of_atoms());

	return 0;
}



