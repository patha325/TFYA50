#ifndef Vec_H
#define Vec_H 

#include <math.h>
#include <iostream>

using namespace std;

class Vec{

private:
	float x;
	float y;
	float z;

public:

	//Constructor
	Vec();
	Vec(float x, float y, float z);

	//Getters
	float getX();
	float getY();
	float getZ();
	float length();

	//Setters
	void setCoords(float x, float y, float z);
	void setX(float x);
	void setY(float y);
	void setZ(float z);

	//Operators
	Vec operator*(float);
	Vec operator+(Vec);
	Vec operator-(Vec);
	Vec& operator+=(Vec&);
	Vec& operator-=(Vec&);
	Vec operator=(Vec);

    //Other functions
    Vec normalize();
    
    //Friends
	friend ostream& operator<<(ostream& os, Vec in);

};

Vec operator*(float,Vec);



#endif