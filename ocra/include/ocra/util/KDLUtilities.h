#ifndef OCRA_UTIL_KDL_UTILITIES_H
#define OCRA_UTIL_KDL_UTILITIES_H

#include <Eigen/Dense>
#include "kdl/frames_io.hpp"
#include "kdl/frames.hpp"
#include <vector>

namespace ocra {
    namespace util {


inline Eigen::VectorXd KDLTwistToEigenVectorXd(const KDL::Twist& kdlTwist)
{
    Eigen::VectorXd tmpVec(6);
    // Linear velocity - Angular Velocity
    tmpVec << kdlTwist(0), kdlTwist(1), kdlTwist(2), kdlTwist(3), kdlTwist(4), kdlTwist(5);
    return tmpVec;
}

inline KDL::Twist EigenVectorToKDLTwist(const Eigen::VectorXd& eigVector)
{
    if (eigVector.size() != 6)
        throw std::runtime_error("[ocra::util::EigenVectorToKDLTwist] wrongly sized Eigen Vector");
    KDL::Vector tmpVel = KDL::Vector(eigVector(0), eigVector(1), eigVector(2));
    KDL::Vector tmpRot = KDL::Vector(eigVector(3), eigVector(4), eigVector(5));
    KDL::Twist tmpTwist = KDL::Twist(tmpVel, tmpRot);
    return tmpTwist;
}
    } // namespace util
} // namespace ocra
#endif // OCRA_UTIL_KDL_UTILITIES_H
