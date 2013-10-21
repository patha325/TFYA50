#include "simulation.h"
#include <iostream>
#include "vec.h"
using namespace std;


int main(){
	//initializegraphics
	cout << "Hej!" << endl;
	float a = 0.2; // Watch out for double warning!
	Simulation* simulation2 = new Simulation(1,1,1,0,0,a,a,a,a,a,a,"Hej!",true);
	simulation2->run_simulation();
	
	system("pause");
	return 0;
}

/*
Origo in the botom left corner, we are in the first octant! all atoms have positive coordinates. (when not moving)





*/