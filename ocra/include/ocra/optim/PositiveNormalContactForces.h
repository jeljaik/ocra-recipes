/** @file PositiveNormalContactForces.h
  * @brief Declaration file of the PositiveNormalContactForces class.
  *
  *
  * @author Jorhabib Eljaik
  *	@date 30/10/2017
  *
  */

#ifndef _POSITIVE_NORMAL_CONTACT_FORCES_H
#define _POSITIVE_NORMAL_CONTACT_FORCES_H

// ocra includes
#include "ocra/optim/LinearFunction.h"

/** @namespace ocra
  * @brief Optimization-based Robot Controller namespace.
  *  a library of classes to write and solve optimization problems dedicated to
  *  the control of multi-body systems.
  */
namespace ocra
{
  /** @class PositiveNormalContactForces
    *	@brief %PositiveNormalContactForces class.
    *	@warning None
    *
    * Implementation of \f$ f_n >= 0\f$.
    */
  class PositiveNormalContactForces : public LinearFunction
  {
    // ------------------------ structures --------------------------------------
  public:
    typedef LinearFunction  functionType_t;     //< alias on the type of the mother class. Needed to duplicate the function tree.

    // ------------------------ constructors ------------------------------------
  private:
    /** Non copyable class */
    //@{
    PositiveNormalContactForces(const PositiveNormalContactForces&);
    PositiveNormalContactForces& operator= (const PositiveNormalContactForces&);
    //@}

  public:
    /** Constructor
      * Builds a linear function Af+b
      */
    PositiveNormalContactForces(Variable& f);

    // ------------------------ public interface --------------------------------
  public:
    /** getters/setters */
    //@{
    //@}

    // ------------------------ private methods ---------------------------------
  private:
    void buildA();
    void buildb();
  };

  // void testPositiveNormalContactForces(void);
}

#endif	//_POSITIVE_NORMAL_CONTACT_FORCES_H

// cmake:sourcegroup=Function
