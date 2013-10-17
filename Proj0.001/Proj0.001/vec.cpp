#include "vec.h"

/*---------------------
CONSTRUCTOR
Saves coordinates x,y,z
---------------------*/
Vec::Vec(float new_x, float new_y, float new_z){

	x = new_x;
	y = new_y;
	z = new_z;
}

/*-----------
GETTERS
-----------*/
float Vec::getX(){
	
	return x;
}
float Vec::getY(){

	return y;
}
float Vec::getZ(){

	return z;
}
float Vec::length(){

	return sqrt(x*x+y*y+z*z);
}

/*----------
SETTERS
----------*/
void Vec::setCoords(float new_x, float new_y, float new_z){

	x = new_x;
	y = new_y;
	z = new_z;
}

void Vec::setX(float new_x){
	
	x = new_x;
}
void Vec::setY(float new_y){

	y = new_y;
}
void Vec::setZ(float new_z){

	z = new_z;
}

/*------------
OPERATORS
------------*/
Vec Vec::operator*(float scalar){

	Vec newVec(this->getX()*scalar,this->getY()*scalar,this->getZ()*scalar);
	return newVec;
}
Vec Vec::operator+(Vec vecA){

	Vec newVec(vecA.getX()+this->getX(),vecA.getY()+this->getY(),vecA.getZ()+this->getZ());
	return newVec;
}

Vec Vec::operator=(Vec vecA){

	x = vecA.getX();
	y = vecA.getY();
	z = vecA.getZ();
	return *this;
}

Vec operator*(float a, Vec vecA){

	return vecA*a;
}