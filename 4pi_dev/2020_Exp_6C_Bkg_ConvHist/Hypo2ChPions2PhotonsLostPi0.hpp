#ifndef _KFCmd_Hypo2ChPions2PhotonsLostPi0_HPP_
#define _KFCmd_Hypo2ChPions2PhotonsLostPi0_HPP_
#include "kfcmd/core/Hypothesis.hpp"


    class Hypo2ChPions2PhotonsLostPi0 : public kfcmd::core::Hypothesis {
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
      Hypo2ChPions2PhotonsLostPi0(double, double, long = 20, double = 1.e-4);
      //! A destructor
      virtual ~Hypo2ChPions2PhotonsLostPi0();
    };

#endif
