// File plumed/src/colvar/Twist2.cpp
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Pbc.h"
#include <string>
#include <cmath>
#include "TwistDerivatives.h"

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR TWIST2
/*
 * TWIST2 colvar is a DNA restraint that keeps a certain twist angle
 * between the chosen base pairs in a DNA fragment.
 * As input provide 48 atoms from bases that are restrained (3*4 bases * 3 atoms)
 * plus 6*2 auxiliary atoms, force const, and desired value of twist.
 * The energy penalty will be added to the potential energy functional
 * E_tw = 0.5 * k * (tw0-tw)^2
 */
//+ENDPLUMEDOC
   
class Twist2 : public Colvar {
  bool pbc;
public:
  explicit Twist2(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Twist2,"TWIST2")

void Twist2::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOMS","Specifies which atoms to use, 6 atoms per base pair");
}

Twist2::Twist2(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if((atoms.size() % 6) != 0)
  {
	  error("Number of specified atoms should be multiple of 6");
	  log.printf("number or atoms %d \n ", atoms.size());
  }

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
  checkRead();
}


// calculator
void Twist2::calculate(){

  if(pbc) makeWhole();

  // Initialize data structure for calculating derivatives.
  Eigen::VectorXd V(getNumberOfAtoms() * 3);  // Flat vector of all coordinates.
  int i = 0;  // Index into flat vector, 3 per atom.
  for(int ia = 0; ia < getNumberOfAtoms(); ia++) {  // ia = index of atom
      const auto position = getPosition(ia);
      V[i++] = position[0];
      V[i++] = position[1];
      V[i++] = position[2];
  }

  // Calculate derivatives.
  auto [twist, gradTwist] = DnaParams::twistDegAndDerivatives(V);

  // Copy derivatives to Plumed.
  i = 0;  // Index into flat vector.
  for(int ia = 0; ia < getNumberOfAtoms(); ia++){
    double x = gradTwist(i++);  
    double y = gradTwist(i++);
    double z = gradTwist(i++);
	  const Vector derivatives(x, y, z);
	  setAtomsDerivatives(ia, derivatives);
   }
      
  setBoxDerivativesNoPbc();
  setValue(twist);

}

}  // namespace
}  // namespace
