#include <TDatabasePDG.h>
#include "Hypo2ChPions2PhotonsLostPi0.hpp"

Hypo2ChPions2PhotonsLostPi0::Hypo2ChPions2PhotonsLostPi0(double energy,
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
  addParticlePxPyPz("pi0", TDatabasePDG::Instance()->GetParticle(111)->Mass());
  addConstantMomentumParticle("origin", energy, Eigen::Vector3d::Zero());
  addEnergyMomentumConstraints("em-vtx0", {getParticle("origin")},
                               {pip, pim, 
				   getParticle("g1"),
				   getParticle("g2"),
				   getParticle("pi0")});
  addOutputVertexConstraintsXYZ("pi+", "vtx0");
  addOutputVertexConstraintsXYZ("pi-", "vtx0");
}

Hypo2ChPions2PhotonsLostPi0::~Hypo2ChPions2PhotonsLostPi0() {}
