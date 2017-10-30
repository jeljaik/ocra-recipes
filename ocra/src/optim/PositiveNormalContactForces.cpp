#include "ocra/optim/PositiveNormalContactForces.h"
#include <math.h>

#include <stdexcept>


namespace ocra
{
  PositiveNormalContactForces::PositiveNormalContactForces(Variable& f)
    :NamedInstance("Positive Normal Contact Forces")
    ,AbilitySet(PARTIAL_X)
    ,CoupledInputOutputSize(false)
    ,LinearFunction(f,1)
  {
    if (f.getSize() != 3)
    {
      std::stringstream ss;
      ss << "[ocra::PositiveNormalContactForces::PositiveNormalContactForces] Size of the variable is not 3 but " << f.getSize();
      throw std::runtime_error(ss.str());
    }

    buildA();
    buildb();
  }

  // ------------------------ private methods ---------------------------------
  void PositiveNormalContactForces::buildA()
  {
     Eigen::Vector3d n;
     n << 0,0,1;
    _jacobian = n.transpose();
  }

  void PositiveNormalContactForces::buildb()
  {
    _b.setConstant(0);
  }

}


//test
namespace ocra
{
  void testPositiveNormalContactForces()
  {
    BaseVariable f("f", 3);
    PositiveNormalContactForces posNormFc(f);

    // Vector3d v;
    //
    // std::cout << coul.getA() << std::endl;
    // std::cout << coul.getb() << std::endl;
    //
    // bool b=true;
    // while (b)
    // {
    //   std::cin >> v[0] >> v[1] >> v[2];
    //   f.setValue(v);
    //   std::cout << coul.getValue() << std::endl;
    //   if (v[0] == -2.71828)
    //     b=false;
    // }
  }
}

// cmake:sourcegroup=Function
