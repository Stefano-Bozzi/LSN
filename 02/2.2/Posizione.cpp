#include "Posizione.h"

//costruttori

//costruttore di default
Posizione::Posizione(){
	m_x=0.;
	m_y=0.;
	m_z=0.;
}

//costruttore a partire da terna cartesiana
Posizione::Posizione(double x, double y, double z){
	m_x=x;
	m_y=y;
	m_z=z;
}

//distruttore di default
Posizione::~Posizione(){
}


// funzioni

//coordinate cartesiane
double Posizione::getX() const {
return m_x;
}
double Posizione::getY() const {
return m_y;
}
double Posizione::getZ() const {
return m_z;
}

//coordinate sferiche
double Posizione::getR() const {
return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
}
double Posizione::getPhi() const {
return atan2(m_y,m_x);
}
double Posizione::getTheta() const {
return acos(m_z/getR());
}

//raggio delle coordinate cilindriche
double Posizione::getRho() const {
	return sqrt(m_x*m_x+m_y*m_y);
}

//distanza da un punto
double Posizione::distx(const Posizione& b) const {
	return sqrt(pow(getX()-b.getX(),2));
}
double Posizione::disty(const Posizione& b) const {
	return sqrt(pow(getY()-b.getY(),2));
}
double Posizione::distz(const Posizione& b) const {
	return sqrt(pow(getZ()-b.getZ(),2));
}
double Posizione::Distanza(const Posizione& b) const {
	return sqrt(pow(getX()-b.getX(),2)+pow(getY()-b.getY(),2)+pow(getZ()-b.getZ(),2));
}
