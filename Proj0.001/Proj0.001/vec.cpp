#include "vec.h"

using namespace std;

/*---------------------
CONSTRUCTOR
Saves coordinates x,y,z
---------------------*/
Vec::Vec(float new_x, float new_y, float new_z){

	x = new_x;
	y = new_y;
	z = new_z;
}

/*-----------------------
CONSTRUCTOR (Alternative)
Constructs a null vector.
------------------------*/
Vec::Vec(){

	x = 0;
	y = 0;
	z = 0;
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
Vec Vec::operator-(Vec vecA){

	Vec newVec(this->getX()-vecA.getX(),this->getY()-vecA.getY(),this->getZ()-vecA.getZ());
	return newVec;
}
Vec& Vec::operator+=(Vec& vecA){

	this->x += vecA.getX();
	this->y += vecA.getY();
	this->z += vecA.getZ();
	return *this;
}
Vec& Vec::operator-=(Vec& vecA){

	this->x -= vecA.getX();
	this->y -= vecA.getY();
	this->z -= vecA.getZ();
	return *this;
}

Vec Vec::operator=(Vec vecA){

	x = vecA.getX();
	y = vecA.getY();
	z = vecA.getZ();
	return *this;
}

bool Vec::operator==(const Vec vecA) const{

    Vec tmp = vecA;
    if (x==tmp.getX() && y==tmp.getY() && z==tmp.getZ()) {
        return true;
    }
    else return false;
}

bool Vec::operator!=(const Vec vecA) const{
    
    Vec tmp = vecA;
    if (x!=tmp.getX() || y!=tmp.getY() || z!=tmp.getZ()) {
        return true;
    }
    else return false;
}

//This function is stupid, I know.
//This is NOT a comparison of vector lengths!
//It needs to look this way for 'maps' to work...
bool Vec::operator<(const Vec vecA) const {
    
    Vec tmp = vecA;
    Vec tmp2 = *this;
    if (tmp2.getX() < tmp.getX()) {
        return true;
    }
    else if((tmp2.getX() == tmp.getX()) && (tmp2.getY() < tmp.getY())){
        return true;
    }
    else if((tmp2.getX() == tmp.getX()) && (tmp2.getY() == tmp.getY()) && (tmp2.getZ() < tmp.getZ())){
        return true;
    }
    else return false;
}

Vec operator*(float a, Vec vecA){

	return vecA*a;
}

/*-------------
 OTHER FUNCTIONS
 -------------*/
Vec Vec::normalize(){
    
    if (length()==0) {
        return Vec(0,0,0);
    }
    else {
        return (1/this->length())*(*this);
    }
}


/*-----------
 FRIENDS
 ----------*/
ostream& operator<<(ostream& os, Vec in){
	os << "(" << in.x << " " << in.y << " "<< in.z << ")";
	return os;
}











