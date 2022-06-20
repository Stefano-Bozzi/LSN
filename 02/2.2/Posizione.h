#ifndef __Posizione_h__
#define __Posizione_h__
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;
class Posizione {

public:

  // costruttori
  Posizione();
  Posizione(double x, double y, double z); 
  // distruttore
  ~Posizione();
  // metodi
  void setX( double x) { m_x = x;};
  void setY( double y) { m_y = y;};  
  void setZ( double z) { m_z = z;};
  double getX() const;       // Coordinate cartesiane
  double getY() const;
  double getZ() const;
  double getR() const;       // Coordinate sferiche
  double getPhi() const;
  double getTheta() const;
  double getRho() const;     // raggio delle coordinate cilindriche
  double distx(const Posizione&) const;
  double disty(const Posizione&) const;
  double distz(const Posizione&) const;
  double Distanza(const Posizione&) const; // distanza da un altro punto

private:

  double m_x, m_y, m_z;  

};

#endif // __posizione_h__
