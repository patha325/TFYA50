using namespace std;
#include "simulation.h"


int main(){
	//initializegraphics

//	Simulation simulation1(0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,"Hej!",true);
	float a = 0.5; // Watch out for double warning!
	Simulation* simulation2 = new Simulation(0,0,0,a,a,a,a,a,a,"Hej!",true);
	simulation2->run_simulation();
	
//	Simulation simulation2(resume.simulation1);


	
	return 0;
}