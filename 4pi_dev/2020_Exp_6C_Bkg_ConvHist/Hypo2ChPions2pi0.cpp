#include <TDatabasePDG.h>
#include "Hypo2ChPions2pi0.hpp"

Hypo2ChPions2pi0::Hypo2ChPions2pi0(double energy,
                                           double magneticField,
                                           long nIter,
                                           double tolerance)
    : kfcmd::core::Hypothesis(energy, magneticField, nIter, tolerance) {
  addVertexXYZ("vtx0");
  auto pip = new kfcmd::core::PiPlusMeson("pi+");
  addChargedParticle(pip);
  auto pim = new kfcmd::core::PiMinusMeson("pi-");
  addChargedParticle(pim);
  addPhoton("g1", "vtx0");
  addPhoton("g2", "vtx0");
  addPhoton("g3", "vtx0");
  addPhoton("g4", "vtx0");
  addConstantMomentumParticle("origin", energy, Eigen::Vector3d::Zero());
  addEnergyMomentumConstraints("em-vtx0", {getParticle("origin")},
                               {pip, pim, 
				   getParticle("g1"),
				   getParticle("g2"),
				   getParticle("g3"),
				   getParticle("g4")});
  addOutputVertexConstraintsXYZ("pi+", "vtx0");
  addOutputVertexConstraintsXYZ("pi-", "vtx0");
  addMassConstraint("pi01Mass", TDatabasePDG::Instance()->GetParticle(111)->Mass(), 
		    {"g1", "g2"});
  addMassConstraint("pi02Mass", TDatabasePDG::Instance()->GetParticle(111)->Mass(), 
		    {"g3", "g4"});
}

Hypo2ChPions2pi0::~Hypo2ChPions2pi0() {}
