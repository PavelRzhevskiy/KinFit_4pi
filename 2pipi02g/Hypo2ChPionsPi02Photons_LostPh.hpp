#ifndef _KFCmd_Hypo2ChPionsPi02Photons_LostPh_HPP_
#define _KFCmd_Hypo2ChPionsPi02Photons_LostPh_HPP_
#include "kfcmd/core/Hypothesis.hpp"

namespace kfcmd {
  namespace hypos {
    /**
     * Implementation of (pi+, pi-, gamma, gamma) hypothesis
     */
    class Hypo2ChPionsPi02Photons_LostPh : public kfcmd::core::Hypothesis {
    public:
      //! A constructor
      /*!
       * @param energy (center-of-mass energy)
       *
       * @param magneticField (magnetic field)
       *
       * @param nIter (maximum number of iterations)
       *
       * @param tolerance (optimization tolerance)
       */
      Hypo2ChPionsPi02Photons_LostPh(double, double, long = 20, double = 1.e-4);
      //! A destructor
      virtual ~Hypo2ChPionsPi02Photons_LostPh();
    };
  } // namespace core
}  // namespace kfcmd

#endif
