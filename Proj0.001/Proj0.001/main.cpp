#include "simulation.h"
#include <iostream>
#include "vec.h"
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
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


int main(){
    
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
	float input_sigma=a;
	float input_epsilon = a;
	float input_lattice_constant = a;
	string input_crystal_structure = "fcc";
	string input_material;
	

	

	
	cout << "Input the number of unit cells in x,y and z direction:" <<endl;
	cin >> input_x;
	cin >> input_y;
	cin >> input_z;
	cout << "Input the wanted time step size:" << endl;
	cin >> input_time_step;
	cout << "Input the wanted number of steps:" << endl;
	cin >> input_steps;
	cout << "Input wanted material (must choose Ar atm):" << endl;
	cin >> input_material;
	get_input_file(input_material,file);
	cout << "Start temperature (K):" << endl;
	cin >> input_temperature;
	cout << "Input cutoff multiples of lattice_constant:" <<endl;
	cin >> input_cutoff;
	input_cutoff = input_cutoff*lattice;

	cout << endl << "------------" << endl;
	cout << "- RUNNING -" << endl;
	cout << "-----------" << endl << endl;
	
	input_sigma = sigma; //[Å]
	input_epsilon = epsilon; //[eV]
	input_crystal_structure = structure;
	input_mass = mass; //[eV/c^2]
	input_lattice_constant=lattice;

	cout << "Sigma: " << sigma << endl;
	cout << "Epsilon: "<<epsilon<<endl;
	cout<<"Structure: "<<structure<<endl;
	cout<<"Mass: "<<mass<<endl;
	cout<<"Lattice: "<<lattice<<endl;


		//input_time_step = 1;
	
	Simulation* simulation2 = new Simulation(input_x,input_y,input_z,input_time_step,input_steps,input_temperature, input_cutoff, 
		input_mass, input_sigma, input_epsilon, input_lattice_constant,input_crystal_structure,false);

	cout << "Running simulation..." << endl << endl;
	simulation2->run_simulation();
	system("pause");
	return 0;
}

/*
Origo in the botom left corner, we are in the first octant! all atoms have positive coordinates. (when not moving)


*/


