#ifndef OCRA_UTIL_KDL_UTILITIES_H
#define OCRA_UTIL_KDL_UTILITIES_H

#include <Eigen/Dense>
#include "kdl/frames_io.hpp"
#include "kdl/frames.hpp"
#include <vector>

namespace ocra {
    namespace util {
        
        /**
         Takes a KDL Twist and converts it into an Eigen Vector. The convention used is the same
         as in KDL, first linear velocity followed by angular velocity.

         @param kdlTwist Input KDL twist.
         @return Corresponding Eigen Vector.
         */
        inline Eigen::VectorXd KDLTwistToEigenVectorXd(const KDL::Twist& kdlTwist)
        {
            Eigen::VectorXd tmpVec(6);
            // Linear velocity - Angular Velocity
            tmpVec << kdlTwist(0), kdlTwist(1), kdlTwist(2), kdlTwist(3), kdlTwist(4), kdlTwist(5);
            return tmpVec;
        }
        
        /**
         Takes a six-dimentional Eigen Vector representing a twist and turns it into a KDL Twist.

         @param eigVector Input Eigen Vector twist.
         @return Corresponding KDL Twist.
         */
        inline KDL::Twist EigenVectorToKDLTwist(const Eigen::VectorXd& eigVector)
        {
            if (eigVector.size() != 6)
                throw std::runtime_error("[ocra::util::EigenVectorToKDLTwist] wrongly sized Eigen Vector");
            KDL::Vector tmpVel = KDL::Vector(eigVector(0), eigVector(1), eigVector(2));
            KDL::Vector tmpRot = KDL::Vector(eigVector(3), eigVector(4), eigVector(5));
            KDL::Twist tmpTwist = KDL::Twist(tmpVel, tmpRot);
            return tmpTwist;
        }

        /**
         Addition between an Eigen Vector and a KDL Twist

         @param eigVector Eigen Vector. Should be of size 6x1.
         @param kdlTwist KDL Twist.
         @return Resulting addition of the two twists.
         */
        template <typename Derived>
        KDL::Twist operator+(const Eigen::DenseBase<Derived> &eigVector, const KDL::Twist &kdlTwist) {
            KDL::Twist kdlTwist1;
            kdlTwist1(0) = kdlTwist(0) + eigVector(0);
            kdlTwist1(1) = kdlTwist(1) + eigVector(1);
            kdlTwist1(2) = kdlTwist(2) + eigVector(2);
            kdlTwist1(3) = kdlTwist(3) + eigVector(3);
            kdlTwist1(4) = kdlTwist(4) + eigVector(4);
            kdlTwist1(5) = kdlTwist(5) + eigVector(5);
            return kdlTwist1;
        }

        /**
         Addition between a KDL Twist and an Eigen Vector

         @param eigVector Eigen Vector. Should be of size 6x1.
         @param kdlTwist KDL Twist.
         @return Resulting addition of the two twists.
         */
        template <typename Derived>
        KDL::Twist operator+(const KDL::Twist &kdlTwist, const Eigen::DenseBase<Derived> &eigVector) {
            KDL::Twist kdlTwist1;
            kdlTwist1(0) = kdlTwist(0) + eigVector(0);
            kdlTwist1(1) = kdlTwist(1) + eigVector(1);
            kdlTwist1(2) = kdlTwist(2) + eigVector(2);
            kdlTwist1(3) = kdlTwist(3) + eigVector(3);
            kdlTwist1(4) = kdlTwist(4) + eigVector(4);
            kdlTwist1(5) = kdlTwist(5) + eigVector(5);
            return kdlTwist1;
        }

        /**
         Multiplication operator between Eigen Matrices and KDL Twists.

         @param eigMatrix Eigen Matrix. Must have 6 columns.
         @param kdlTwist KDL Twist.
         @return Result of the multiplication as a KDL Twist.
         */
        template <typename Derived>
        KDL::Twist operator*(const Eigen::MatrixBase<Derived> &eigMatrix, const KDL::Twist &kdlTwist) {
            Eigen::VectorXd tmpEig(6);
            tmpEig(0) = kdlTwist(0);
            tmpEig(1) = kdlTwist(1);
            tmpEig(2) = kdlTwist(2);
            tmpEig(3) = kdlTwist(3);
            tmpEig(4) = kdlTwist(4);
            tmpEig(5) = kdlTwist(5);

            Eigen::VectorXd res(6);
            res = eigMatrix*tmpEig;
            KDL::Twist outPut = KDL::Twist(KDL::Vector(res(0), res(1), res(2)), KDL::Vector(res(3), res(4), res(5)));
            return outPut;
        }
        
        /**
         Multiplication operator between Eigen Matrices and KDL Twists.
         
         @param eigMatrix Eigen Matrix. Must have 6 columns.
         @param kdlTwist KDL Twist.
         @return Result of the multiplication as a KDL Twist.
         */
        template <typename Derived>
        KDL::Twist operator*(Eigen::MatrixBase<Derived> &eigMatrix, const KDL::Twist &kdlTwist)  {
            Eigen::VectorXd tmpEig(6);
            tmpEig(0) = kdlTwist(0);
            tmpEig(1) = kdlTwist(1);
            tmpEig(2) = kdlTwist(2);
            tmpEig(3) = kdlTwist(3);
            tmpEig(4) = kdlTwist(4);
            tmpEig(5) = kdlTwist(5);
            
            Eigen::VectorXd res(6);
            res = eigMatrix*tmpEig;
            KDL::Twist outPut = KDL::Twist(KDL::Vector(res(0), res(1), res(2)), KDL::Vector(res(3), res(4), res(5)));
            return outPut;
        }
        
        /**
         Computes the skew symmetric matrix from an input KDL frame.
         
         @param inputFrame KDL Frame.
         @return Skew symmetric matrix.
         */
        inline Eigen::Matrix3d computeSkewSymmetric(const KDL::Frame &inputFrame) {
            // Allocate space for output skew matrix
            Eigen::Matrix3d skew;
            // Extract position vector from kDL frame
            Eigen::Vector3d p; p << inputFrame.p.x(), inputFrame.p.y(), inputFrame.p.z();
            skew <<    0,  -p(2),   p(1),
                    p(2),      0,  -p(0),
                   -p(1),   p(0),      0;
            return skew;
        }

        
        /**
         Takes a KDL Frame and transforms it into an Eigen 4x4 matrix.
         
         @param  input Input KDL rotation matrix.
         @return output Output Eigen 3D matrix.
         */
        inline Eigen::MatrixXd KDLFrameToEigenHomogeneous(const KDL::Frame &input) {
            Eigen::MatrixXd output(4,4);
            KDL::Frame copyInput = input;
            Eigen::VectorXd H_vec(16);
            copyInput.Make4x4(H_vec.data());
            output << H_vec;
            output.transposeInPlace();
            return output;
        }
        
        /**
         Gets the adjoint of a KDL::Frame, i.e. given
           adj = [R       0;
                 S(p)*R   R]
         
         for H = [R  p;
                  0  1]
         
          where S(P) is the skew symmetric matrix corresponding to the cross
          product operation.

         @param[in] inputFrame Input KDL Frame whose blocks will be used to compute the adjoint transformation.
         @return    Adjoint transform from the input frame.
         */
        inline Eigen::MatrixXd getKDLFrameAdjoint(const KDL::Frame &inputFrame) {
            // Allocate space for adjoint
            Eigen::MatrixXd adjoint(6,6);
            adjoint.setZero();
            // Extract homogeneous matrix from KDL Frame
            Eigen::MatrixXd tmpHomogeneous(4,4);
            tmpHomogeneous = KDLFrameToEigenHomogeneous(inputFrame);
            // Build the adjoint matrix
            // First copy the rotation block
            adjoint.block(0,0,3,3) = tmpHomogeneous.block(0,0,3,3);
            adjoint.block(3,3,3,3) = tmpHomogeneous.block(0,0,3,3);
            adjoint.block(3,0,3,3) = computeSkewSymmetric(inputFrame) * tmpHomogeneous.block(0,0,3,3);
            return adjoint;
            
        }
        
        // /**
        //  Returns pointer to data stored in KDL Twist. Well, actually to the copy of this data.
        //
        //  @param kdlTwist Input KDL Twist.
        //  @return Pointer to data.
        //  */
        // double* getTwistData(const KDL::Twist &kdlTwist) {
        //     std::vector<double> outTwist(3);
        //     outTwist.assign(&kdlTwist.vel.data[0], &kdlTwist.vel.data[0] + 3);
        //     outTwist.insert(outTwist.end(), &kdlTwist.rot.data[0], &kdlTwist.rot.data[0] + 3);
        //     return outTwist.data();
        // }

    } // namespace util
} // namespace ocra
#endif // OCRA_UTIL_KDL_UTILITIES_H
