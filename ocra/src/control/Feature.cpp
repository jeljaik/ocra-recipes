#include "ocra/control/Feature.h"


#include <iostream>

// TODO: Clean up old code to convert frame of references to the link frame rather than in the world frame. Note: all frames should be expressed in the world frame.

using namespace Eigen;

namespace ocra
{
  // --- ABSTRACT -----------------------------------------------

  Feature::Feature(const std::string& name)
    : NamedInstance(name)
  {
  }

  Feature::~Feature()
  {
  }


  // --- POSITION -----------------------------------------------

  struct PositionFeature::Pimpl
  {
    ControlFrame::Ptr controlFrame;
    ECartesianDof axes;
    MatrixXd u;
    MatrixXd spaceTransform;
    VectorXd error;
    VectorXd errorDot;
    VectorXd acceleration;
    VectorXd effort;
    MatrixXd jacobian;
    MatrixXd M;
    MatrixXd Minv;

    Pimpl(ControlFrame::Ptr cf, ECartesianDof a)
      : controlFrame(cf), axes(a)
    {
      int dim = utils::computeDimensionFor(axes, NONE);

      u.resize(3, dim);
      int k = 0;
      if(axes & X)
        u.col(k++) = Vector3d(1., 0., 0.);
      if(axes & Y)
        u.col(k++) = Vector3d(0., 1., 0.);
      if(axes & Z)
        u.col(k++) = Vector3d(0., 0., 1.);

      error = VectorXd::Zero(dim);
      errorDot = VectorXd::Zero(dim);
      effort = VectorXd::Zero(dim);
      acceleration = VectorXd::Zero(dim);
      jacobian = MatrixXd::Zero(dim, cf->getJacobian().cols());
    }
  };

  PositionFeature::PositionFeature(const std::string& name, ControlFrame::Ptr frame, ECartesianDof axes)
    : Feature(name)
    , pimpl(new Pimpl(frame, axes))
  {

  }

  const MatrixXd& PositionFeature::getSpaceTransform() const
  {
#ifndef OCRA_USES_KDL
    const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
    pimpl->spaceTransform = pimpl->u.transpose() * R.adjoint();
#else
    //TODO: Review the assumption that the adjoint of a quaternion in LGSM is equivalent to the rotation matrix of a KDL frame.
    Eigen::Matrix3d controlFrameRotation = ocra::util::rotationMatrixFromKDLFrame(pimpl->controlFrame->getPositionKDL());
    pimpl->spaceTransform = pimpl->u.transpose() * controlFrameRotation;
#endif
    return pimpl->spaceTransform;
  }

  int PositionFeature::getDimension() const
  {
    return utils::computeDimensionFor(pimpl->axes, NONE);
  }

  const VectorXd& PositionFeature::computeEffort(const Feature& featureDes) const
  {
    const PositionFeature& sdes = dynamic_cast<const PositionFeature&>(featureDes);
#ifndef OCRA_USES_KDL

//    // first compute the linear velocity error in the mobile frame
//    const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
//    const Eigen::Displacementd::Rotation3D& Rdes = sdes.pimpl->controlFrame->getPosition().getRotation();
//    const Eigen::Displacementd::Rotation3D Rdes_in_r = R.inverse() * Rdes;
//    const Vector3d eff = pimpl->controlFrame->getWrench().getForce() - Rdes_in_r.adjoint() * sdes.pimpl->controlFrame->getWrench().getForce();

//    // then project it on the controlled axes
//    pimpl->effort = getSpaceTransform() * eff;

    pimpl->effort = pimpl->u.transpose() * (pimpl->controlFrame->getWrench().getForce() - sdes.pimpl->controlFrame->getWrench().getForce());
#else
    KDL::Vector forceDiff = pimpl->controlFrame->getWrenchKDL().force - sdes.pimpl->controlFrame->getWrenchKDL().force;
    pimpl->effort = pimpl->u.transpose() * ocra::util::KDLVectorToEigenVector3d(forceDiff);
#endif
    return pimpl->effort;
  }

  const VectorXd& PositionFeature::computeEffort() const
  {
#ifndef OCRA_USES_KDL
//    // first compute the linear velocity error in the mobile frame
//    const Vector3d eff = pimpl->controlFrame->getWrench().getForce();

//    // then project it on the controlled axes
//    pimpl->effort = getSpaceTransform() * eff;

    pimpl->effort = pimpl->u.transpose() * pimpl->controlFrame->getWrench().getForce();
#else
    pimpl->effort = pimpl->u.transpose() * ocra::util::KDLVectorToEigenVector3d(pimpl->controlFrame->getWrenchKDL().force);
#endif
    return pimpl->effort;
  }

  const VectorXd& PositionFeature::computeAcceleration(const Feature& featureDes) const
  {
    const PositionFeature& sdes = dynamic_cast<const PositionFeature&>(featureDes);
#ifndef OCRA_USES_KDL
//    const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
//    const Eigen::Displacementd::Rotation3D& Rdes = sdes.pimpl->controlFrame->getPosition().getRotation();
//    const Eigen::Displacementd::Rotation3D Rdes_in_r = R.inverse() * Rdes;
//    const VectorXd acc = pimpl->controlFrame->getAcceleration().getLinearVelocity() - Rdes_in_r.adjoint() * sdes.pimpl->controlFrame->getAcceleration().getLinearVelocity();

//    pimpl->acceleration = getSpaceTransform() * acc;

    pimpl->acceleration = pimpl->u.transpose() * (pimpl->controlFrame->getAcceleration().getLinearVelocity() - sdes.pimpl->controlFrame->getAcceleration().getLinearVelocity());
#else
    Eigen::Vector3d accDiff = ocra::util::KDLVectorToEigenVector3d(pimpl->controlFrame->getAccelerationKDL().vel - sdes.pimpl->controlFrame->getAccelerationKDL().vel);
    pimpl->acceleration = pimpl->u.transpose() * (accDiff);
#endif
    return pimpl->acceleration;
  }

  const VectorXd& PositionFeature::computeAcceleration() const
  {
#ifndef OCRA_USES_KDL
//    const VectorXd acc = pimpl->controlFrame->getAcceleration().getLinearVelocity();
//    pimpl->acceleration = getSpaceTransform() * acc;

    pimpl->acceleration = pimpl->u.transpose() * pimpl->controlFrame->getAcceleration().getLinearVelocity();
#else
    Eigen::Vector3d linAcc = ocra::util::KDLVectorToEigenVector3d(pimpl->controlFrame->getAccelerationKDL().vel);
    pimpl->acceleration = pimpl->u.transpose() * linAcc;
#endif
    return pimpl->acceleration;
  }

  const VectorXd& PositionFeature::computeError(const Feature& featureDes) const
  {
    const PositionFeature& sdes = dynamic_cast<const PositionFeature&>(featureDes);
#ifndef OCRA_USES_KDL
////    // first compute the position error in the mobile frame
//    const Vector3d e0 = pimpl->controlFrame->getPosition().getTranslation() - sdes.pimpl->controlFrame->getPosition().getTranslation();
//    const Vector3d e_in_r = pimpl->controlFrame->getPosition().getRotation().inverse() * e0;

//    // then project it on the controlled axes
//    pimpl->error = getSpaceTransform() * e_in_r;

    pimpl->error = pimpl->u.transpose() * (pimpl->controlFrame->getPosition().getTranslation() - sdes.pimpl->controlFrame->getPosition().getTranslation());
#else
    KDL::Vector transDiffKDL = pimpl->controlFrame->getPositionKDL().p - sdes.pimpl->controlFrame->getPositionKDL().p;
    Eigen::Vector3d transDiff = ocra::util::KDLVectorToEigenVector3d(transDiffKDL);
    pimpl->error = pimpl->u.transpose() * transDiff;
#endif
    return pimpl->error;
  }

  const VectorXd& PositionFeature::computeError() const
  {
#ifndef OCRA_USES_KDL
//    // first compute the position error in the mobile frame
//    const Vector3d e0 = pimpl->controlFrame->getPosition().getTranslation();
//    const Vector3d e_in_r = pimpl->controlFrame->getPosition().getRotation().inverse() * e0;

//    // then project it on the controlled axes
//    pimpl->error = getSpaceTransform() * e_in_r;

    pimpl->error = pimpl->u.transpose() * pimpl->controlFrame->getPosition().getTranslation();
#else
    Eigen::Vector3d tmp = ocra::util::KDLVectorToEigenVector3d(pimpl->controlFrame->getPositionKDL().p);
    pimpl->error = pimpl->u.transpose() * tmp;
#endif
    return pimpl->error;
  }

  const VectorXd& PositionFeature::computeErrorDot(const Feature& featureDes) const
  {
    const PositionFeature& sdes = dynamic_cast<const PositionFeature&>(featureDes);
#ifndef OCRA_USES_KDL
//    // first compute the linear velocity error in the mobile frame
//    const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
//    const Eigen::Displacementd::Rotation3D& Rdes = sdes.pimpl->controlFrame->getPosition().getRotation();
//    const Eigen::Displacementd::Rotation3D Rdes_in_r = R.inverse() * Rdes;
//    const Vector3d errDot = pimpl->controlFrame->getVelocity().getLinearVelocity() - Rdes_in_r.adjoint() * sdes.pimpl->controlFrame->getVelocity().getLinearVelocity();


//    // then project it on the controlled axes
//    pimpl->errorDot = getSpaceTransform() * errDot;

    pimpl->errorDot = pimpl->u.transpose() * (pimpl->controlFrame->getVelocity().getLinearVelocity()-sdes.pimpl->controlFrame->getVelocity().getLinearVelocity());
#else
    KDL::Vector linVelKDL = pimpl->controlFrame->getVelocityKDL().vel - sdes.pimpl->controlFrame->getVelocityKDL().vel;
    Eigen::Vector3d linVel = ocra::util::KDLVectorToEigenVector3d(linVelKDL);
    pimpl->errorDot = pimpl->u.transpose() * linVel;
#endif
    return pimpl->errorDot;
  }

  const VectorXd& PositionFeature::computeErrorDot() const
  {
#ifndef OCRA_USES_KDL
//    // first compute the linear velocity error in the mobile frame
//    const Vector3d errDot = pimpl->controlFrame->getVelocity().getLinearVelocity();

//    // then project it on the controlled axes
//    pimpl->errorDot = getSpaceTransform() * errDot;

    pimpl->errorDot = pimpl->u.transpose() * pimpl->controlFrame->getVelocity().getLinearVelocity();

#else
    KDL::Vector tmpLinVel = pimpl->controlFrame->getVelocityKDL().vel;
    pimpl->errorDot = pimpl->u.transpose() * ocra::util::KDLVectorToEigenVector3d(tmpLinVel);
#endif
    return pimpl->errorDot;
  }

  const MatrixXd& PositionFeature::computeJacobian(const Feature& featureDes) const
  {
//    const PositionFeature& sdes = dynamic_cast<const PositionFeature&>(featureDes);

//    // first compute the jacobian of the position error in the mobile frame
//    const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
//    const Eigen::Displacementd::Rotation3D& Rdes = sdes.pimpl->controlFrame->getPosition().getRotation();
//    const Eigen::Displacementd::Rotation3D Rdes_in_r = R.inverse() * Rdes;
//    const MatrixXd jacobian = pimpl->controlFrame->getJacobian().bottomRows(3) - Rdes_in_r.adjoint() * sdes.pimpl->controlFrame->getJacobian().bottomRows(3);


    // then project it on the controlled axes
//    pimpl->jacobian = getSpaceTransform() * jacobian;

    pimpl->jacobian = pimpl->u.transpose() * pimpl->controlFrame->getJacobian().bottomRows(3);

    return pimpl->jacobian;
  }

  const MatrixXd& PositionFeature::computeJacobian() const
  {
//    // first compute the jacobian of the position error in the mobile frame
//    const MatrixXd jacobian = pimpl->controlFrame->getJacobian().bottomRows(3);

//    // then project it on the controlled axes
//    pimpl->jacobian = getSpaceTransform() * jacobian;

    pimpl->jacobian = pimpl->u.transpose() * pimpl->controlFrame->getJacobian().bottomRows(3);
    return pimpl->jacobian;
  }

  const Eigen::MatrixXd& PositionFeature::computeProjectedMass(const Feature& featureDes) const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[PositionFeature::computeProjectedMass] feature does not depend on configuration");

    pimpl->M = computeProjectedMassInverse(featureDes).inverse();

    return pimpl->M;
  }

  const Eigen::MatrixXd& PositionFeature::computeProjectedMass() const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[PositionFeature::computeProjectedMass] feature does not depend on configuration");

    pimpl->M = computeProjectedMassInverse().inverse();

    return pimpl->M;
  }

  const Eigen::MatrixXd& PositionFeature::computeProjectedMassInverse(const Feature& featureDes) const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[PositionFeature::computeProjectedMassInverse] feature does not depend on configuration");

    const MatrixXd& J = computeJacobian(featureDes);
    const MatrixXd& Minv = pimpl->controlFrame->getModel().getInertiaMatrixInverse();
    pimpl->Minv = J * Minv * J.transpose();

    return pimpl->Minv;
  }

  const Eigen::MatrixXd& PositionFeature::computeProjectedMassInverse() const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[PositionFeature::computeProjectedMassInverse] feature does not depend on configuration");

    const MatrixXd& J = computeJacobian();
    const MatrixXd& Minv = pimpl->controlFrame->getModel().getInertiaMatrixInverse();
    pimpl->Minv = J * Minv * J.transpose();

    return pimpl->Minv;
  }

  TaskState PositionFeature::getState() const
  {
    TaskState state;
#ifndef OCRA_USES_KDL
    state.setPosition(pimpl->controlFrame->getPosition());
    state.setVelocity(pimpl->controlFrame->getVelocity());
    //   state.setAcceleration(pimpl->controlFrame->getAcceleration());
    //   state.setWrench(pimpl->controlFrame->getWrench());
#else
    state.setPositionKDL(pimpl->controlFrame->getPositionKDL());
    state.setVelocityKDL(pimpl->controlFrame->getVelocityKDL());
#endif
      return state;
  }

  void PositionFeature::setState(const TaskState& newState)
  {
#ifndef OCRA_USES_KDL
      try {
          TargetFrame::Ptr targetFrame = std::dynamic_pointer_cast<TargetFrame>(pimpl->controlFrame);
          if(newState.hasPosition()) {
                targetFrame->setPosition(newState.getPosition());
            }
            if(newState.hasVelocity()) {
                targetFrame->setVelocity(newState.getVelocity());
            }
            if(newState.hasAcceleration()) {
                targetFrame->setAcceleration(newState.getAcceleration());
            }
            if(newState.hasWrench()) {
                targetFrame->setWrench(newState.getWrench());
            }
      } catch (int errCode) {
          std::cout << "You cannot set the state of this feature because it is not a desired feature. It must be constructed with a TargetFrame." << errCode << std::endl;
      }
#else
      try {
          TargetFrame::Ptr targetFrame = std::dynamic_pointer_cast<TargetFrame>(pimpl->controlFrame);
          if(newState.hasPositionKDL()) {
              targetFrame->setPositionKDL(newState.getPositionKDL());
          }
          if(newState.hasVelocityKDL()) {
              targetFrame->setVelocityKDL(newState.getVelocityKDL());
          }
          if(newState.hasAccelerationKDL()) {
              targetFrame->setAccelerationKDL(newState.getAccelerationKDL());
          }
          if(newState.hasWrenchKDL()) {
              targetFrame->setWrenchKDL(newState.getWrenchKDL());
          }
      } catch (int errCode) {
          std::cout << "You cannot set the state of this feature because it is not a desired feature. It must be constructed with a TargetFrame." << errCode << std::endl;
      }
#endif
  }


  // --- POINT CONTACT ------------------------------------------

  struct PointContactFeature::Pimpl
  {
    ControlFrame::Ptr controlFrame;
    MatrixXd spaceTransform;
    VectorXd error;
    VectorXd errorDot;
    VectorXd effort;
    VectorXd acceleration;
    MatrixXd jacobian;
    MatrixXd M;
    MatrixXd Minv;

    Pimpl(ControlFrame::Ptr cf)
      : controlFrame(cf)
      , spaceTransform(Matrix3d::Identity())
      , error(VectorXd::Zero(3))
      , errorDot(VectorXd::Zero(3))
      , effort(VectorXd::Zero(3))
      , acceleration(VectorXd::Zero(3))
      , jacobian(MatrixXd::Zero(3, cf->getJacobian().cols()))
    {
    }
  };

  PointContactFeature::PointContactFeature(const std::string& name, ControlFrame::Ptr frame)
    : Feature(name)
    , pimpl(new Pimpl(frame))
  {
  }

  const MatrixXd& PointContactFeature::getSpaceTransform() const
  {
    return pimpl->spaceTransform;
  }

  int PointContactFeature::getDimension() const
  {
    return 3;
  }

  const VectorXd& PointContactFeature::computeEffort(const Feature& featureDes) const
  {
    throw std::runtime_error("[PointContactFeature::computeEffort(const Feature&)] Desired feature are irrelevant in PointContactFeatures");
  }

  const VectorXd& PointContactFeature::computeEffort() const
  {
    pimpl->effort = pimpl->controlFrame->getWrench().getForce();
    return pimpl->effort;
  }

  const VectorXd& PointContactFeature::computeAcceleration(const Feature& featureDes) const
  {
    throw std::runtime_error("[PointContactFeature::computeAcceleration(const Feature&)] Desired feature are irrelevant in PointContactFeatures");
  }

  const VectorXd& PointContactFeature::computeAcceleration() const
  {
    pimpl->acceleration = pimpl->controlFrame->getAcceleration().getLinearVelocity();
    return pimpl->acceleration;
  }

  const VectorXd& PointContactFeature::computeError(const Feature& featureDes) const
  {
    throw std::runtime_error("[PointContactFeature::computeError(const Feature&)] Desired feature are irrelevant in PointContactFeatures");
  }

  const VectorXd& PointContactFeature::computeError() const
  {
//    const Vector3d& e0 = pimpl->controlFrame->getPosition().getTranslation();
//    pimpl->error = pimpl->controlFrame->getPosition().getRotation().inverse() * e0;

    pimpl->error = pimpl->controlFrame->getPosition().getTranslation();
    return pimpl->error;
  }

  const VectorXd& PointContactFeature::computeErrorDot(const Feature& featureDes) const
  {
    throw std::runtime_error("[PointContactFeature::computeErrorDot(const Feature&)] Desired feature are irrelevant in PointContactFeatures");
  }

  const VectorXd& PointContactFeature::computeErrorDot() const
  {
    pimpl->errorDot = pimpl->controlFrame->getVelocity().getLinearVelocity();
    return pimpl->errorDot;
  }

  const MatrixXd& PointContactFeature::computeJacobian(const Feature& featureDes) const
  {
    throw std::runtime_error("[PointContactFeature::computeJacobian(const Feature&)] Desired feature are irrelevant in PointContactFeatures");
  }

  const MatrixXd& PointContactFeature::computeJacobian() const
  {
    pimpl->jacobian = pimpl->controlFrame->getJacobian().bottomRows(3);
    return pimpl->jacobian;
  }

  const Eigen::MatrixXd& PointContactFeature::computeProjectedMass(const Feature& featureDes) const
  {
    throw std::runtime_error("[PointContactFeature::computeProjectedMass(const Feature&)] Desired feature are irrelevant in PointContactFeatures");
  }

  const Eigen::MatrixXd& PointContactFeature::computeProjectedMass() const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[PositionFeature::computeProjectedMass] Feature must depend on configuration!");

    pimpl->M = computeProjectedMassInverse().inverse();

    return pimpl->M;
  }

  const Eigen::MatrixXd& PointContactFeature::computeProjectedMassInverse(const Feature& featureDes) const
  {
    throw std::runtime_error("[PointContactFeature::computeProjectedMassInverse(const Feature&)] Desired feature are irrelevant in PointContactFeatures");
  }

  const Eigen::MatrixXd& PointContactFeature::computeProjectedMassInverse() const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[PositionFeature::computeProjectedMassInverse] feature does not depend on configuration");

    const MatrixXd& J = computeJacobian();
    const MatrixXd& Minv = pimpl->controlFrame->getModel().getInertiaMatrixInverse();
    pimpl->Minv = J * Minv * J.transpose();

    return pimpl->Minv;
  }
  TaskState PointContactFeature::getState() const
  {
      TaskState state;
      state.setPosition(pimpl->controlFrame->getPosition());
      state.setVelocity(pimpl->controlFrame->getVelocity());
      state.setAcceleration(pimpl->controlFrame->getAcceleration());
      state.setWrench(pimpl->controlFrame->getWrench());

      return state;
  }

  void PointContactFeature::setState(const TaskState& newState)
  {
      try {
          TargetFrame::Ptr targetFrame = std::dynamic_pointer_cast<TargetFrame>(pimpl->controlFrame);
          if(newState.hasPosition()) {
                targetFrame->setPosition(newState.getPosition());
            }
            if(newState.hasVelocity()) {
                targetFrame->setVelocity(newState.getVelocity());
            }
            if(newState.hasAcceleration()) {
                targetFrame->setAcceleration(newState.getAcceleration());
            }
            if(newState.hasWrench()) {
                targetFrame->setWrench(newState.getWrench());
            }
      } catch (int errCode) {
          std::cout << "You cannot set the state of this feature because it is not a desired feature. It must be constructed with a TargetFrame." << errCode << std::endl;
      }
  }

  // --- ORIENTATION --------------------------------------------

  struct OrientationFeature::Pimpl
  {
    ControlFrame::Ptr controlFrame;
    MatrixXd spaceTransform;
    VectorXd error;
    VectorXd errorDot;
    VectorXd effort;
    VectorXd acceleration;
    MatrixXd jacobian;
    MatrixXd M;
    MatrixXd Minv;

    Pimpl(ControlFrame::Ptr cf)
      : controlFrame(cf)
    {
      spaceTransform = MatrixXd::Identity(3, 3);
      error = VectorXd::Zero(3);
      errorDot = VectorXd::Zero(3);
      effort = VectorXd::Zero(3);
      acceleration = VectorXd::Zero(3);
      jacobian = MatrixXd::Zero(3, cf->getJacobian().cols());
    }
  };

  OrientationFeature::OrientationFeature(const std::string& name, ControlFrame::Ptr frame)
    : Feature(name)
    , pimpl(new Pimpl(frame))
  {
  }

  const MatrixXd& OrientationFeature::getSpaceTransform() const
  {
    return pimpl->spaceTransform;
  }

  int OrientationFeature::getDimension() const
  {
    return 3;
  }

  const VectorXd& OrientationFeature::computeEffort(const Feature& featureDes) const
  {
    const OrientationFeature& sdes = dynamic_cast<const OrientationFeature&>(featureDes);
#ifndef OCRA_USES_KDL

    const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
    const Eigen::Displacementd::Rotation3D& Rdes = sdes.pimpl->controlFrame->getPosition().getRotation();
    const Eigen::Displacementd::Rotation3D Rdes_in_r = R.inverse() * Rdes;

    pimpl->effort = pimpl->controlFrame->getWrench().getTorque() - Rdes_in_r.adjoint() * sdes.pimpl->controlFrame->getWrench().getTorque();
#else
    using namespace ocra::util;

    const Eigen::Matrix3d& R = rotationMatrixFromKDLFrame(pimpl->controlFrame->getPositionKDL());
    const Eigen::Matrix3d& Rdes = rotationMatrixFromKDLFrame(sdes.pimpl->controlFrame->getPositionKDL());
    const Eigen::Matrix3d Rdes_in_r = R.inverse() * Rdes;

    pimpl->effort = KDLVectorToEigenVector3d(pimpl->controlFrame->getWrenchKDL().torque) - Rdes_in_r * KDLVectorToEigenVector3d(sdes.pimpl->controlFrame->getWrenchKDL().torque);
#endif
    return pimpl->effort;
  }

  const VectorXd& OrientationFeature::computeEffort() const
  {
#ifndef OCRA_USES_KDL
    pimpl->effort = pimpl->controlFrame->getWrench().getTorque();
#else
    pimpl->effort = ocra::util::KDLVectorToEigenVector3d(pimpl->controlFrame->getWrenchKDL().torque);
#endif
    return pimpl->effort;
  }

  const VectorXd& OrientationFeature::computeAcceleration(const Feature& featureDes) const
  {
    const OrientationFeature& sdes = dynamic_cast<const OrientationFeature&>(featureDes);
#ifndef OCRA_USES_KDL
    const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
//    const Eigen::Displacementd::Rotation3D& Rdes = sdes.pimpl->controlFrame->getPosition().getRotation();
//    const Eigen::Displacementd::Rotation3D Rdes_in_r = R.inverse() * Rdes;

//    pimpl->acceleration = pimpl->controlFrame->getAcceleration().getAngularVelocity() - Rdes_in_r.adjoint() * sdes.pimpl->controlFrame->getAcceleration().getAngularVelocity();

    pimpl->acceleration = pimpl->controlFrame->getAcceleration().getAngularVelocity() - sdes.pimpl->controlFrame->getAcceleration().getAngularVelocity();
#else
    using namespace ocra::util;
    pimpl->acceleration = KDLVectorToEigenVector3d(pimpl->controlFrame->getAccelerationKDL().rot - sdes.pimpl->controlFrame->getAccelerationKDL().rot);
#endif
    return pimpl->acceleration;
  }

  const VectorXd& OrientationFeature::computeAcceleration() const
  {
#ifndef OCRA_USES_KDL
    pimpl->acceleration = pimpl->controlFrame->getAcceleration().getAngularVelocity();
#else
    pimpl->acceleration = ocra::util::KDLVectorToEigenVector3d(pimpl->controlFrame->getAccelerationKDL().rot);
#endif
    return pimpl->acceleration;
  }

  const VectorXd& OrientationFeature::computeError(const Feature& featureDes) const
  {
    const OrientationFeature& sdes = dynamic_cast<const OrientationFeature&>(featureDes);
#ifndef OCRA_USES_KDL
    const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
    const Eigen::Displacementd::Rotation3D& Rdes = sdes.pimpl->controlFrame->getPosition().getRotation();

//    pimpl->error = (Rdes.inverse() * R).log();
    pimpl->error = Rdes.adjoint()*((Rdes.inverse() * R).log());
#else
    //TODO: REVIEW! wr.t. OrientationFeature::computeError.
    using namespace ocra::util;
    //        const Eigen::Matrix3d& R = rotationMatrixFromKDLFrame(pimpl->controlFrame->getPositionKDL());
    const Eigen::Matrix3d& Rdes = rotationMatrixFromKDLFrame(sdes.pimpl->controlFrame->getPositionKDL());

    const KDL::Frame& frame = pimpl->controlFrame->getPositionKDL();
    const KDL::Frame& frameDes = sdes.pimpl->controlFrame->getPositionKDL();
    const KDL::Frame tmp = frameDes.Inverse()*frame;
    Eigen::Vector3d tmpLog;
    quaternionLogFromKDLFrame(tmp, tmpLog);
    pimpl->error = Rdes*tmpLog;
#endif
    return pimpl->error;
  }

  const VectorXd& OrientationFeature::computeError() const
  {
#ifndef OCRA_USES_KDL
    const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
    pimpl->error = R.log();
#else
    //TODO: REVIEW!! w.r.t OrientationFeature::computeError()
    const KDL::Frame& tmp = pimpl->controlFrame->getPositionKDL();
    Eigen::Vector3d tmpLog;
    ocra::util::quaternionLogFromKDLFrame(tmp, tmpLog);
    pimpl->error = tmpLog;
#endif
    return pimpl->error;
  }

  const VectorXd& OrientationFeature::computeErrorDot(const Feature& featureDes) const
  {
    const OrientationFeature& sdes = dynamic_cast<const OrientationFeature&>(featureDes);
#ifndef OCRA_USES_KDL
    const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
    const Eigen::Displacementd::Rotation3D& Rdes = sdes.pimpl->controlFrame->getPosition().getRotation();
    const Eigen::Displacementd::Rotation3D Rdes_in_r = R.inverse() * Rdes;

//    pimpl->errorDot = pimpl->controlFrame->getVelocity().getAngularVelocity() - Rdes_in_r.adjoint() * sdes.pimpl->controlFrame->getVelocity().getAngularVelocity();
    pimpl->errorDot = pimpl->controlFrame->getVelocity().getAngularVelocity() - sdes.pimpl->controlFrame->getVelocity().getAngularVelocity();
#else
    KDL::Vector velDiff = pimpl->controlFrame->getVelocityKDL().rot - sdes.pimpl->controlFrame->getVelocityKDL().rot;
    pimpl->errorDot = ocra::util::KDLVectorToEigenVector3d(velDiff);
#endif
    return pimpl->errorDot;
  }

  const VectorXd& OrientationFeature::computeErrorDot() const
  {
#ifndef OCRA_USES_KDL
    pimpl->errorDot = pimpl->controlFrame->getVelocity().getAngularVelocity();
#else
    pimpl->errorDot = ocra::util::KDLVectorToEigenVector3d(pimpl->controlFrame->getVelocityKDL().rot);
#endif
    return pimpl->errorDot;
  }

  const MatrixXd& OrientationFeature::computeJacobian(const Feature& featureDes) const
  {
    const OrientationFeature& sdes = dynamic_cast<const OrientationFeature&>(featureDes);

    const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
    const Eigen::Displacementd::Rotation3D& Rdes = sdes.pimpl->controlFrame->getPosition().getRotation();
    const Eigen::Displacementd::Rotation3D Rdes_in_r = R.inverse() * Rdes;

//    pimpl->jacobian = pimpl->controlFrame->getJacobian().topRows(3) - Rdes_in_r.adjoint() * sdes.pimpl->controlFrame->getJacobian().topRows(3);

    pimpl->jacobian = pimpl->controlFrame->getJacobian().topRows(3) - sdes.pimpl->controlFrame->getJacobian().topRows(3);
    return pimpl->jacobian;
  }

  const MatrixXd& OrientationFeature::computeJacobian() const
  {
    pimpl->jacobian = pimpl->controlFrame->getJacobian().topRows(3);
    return pimpl->jacobian;
  }

  const MatrixXd& OrientationFeature::computeProjectedMass(const Feature& featureDes) const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[OrientationFeature::computeProjectedMass] feature does not depend on configuration");

    pimpl->M = computeProjectedMassInverse(featureDes).inverse();

    return pimpl->M;
  }

  const MatrixXd& OrientationFeature::computeProjectedMass() const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[OrientationFeature::computeProjectedMass] feature does not depend on configuration");

    pimpl->M = computeProjectedMassInverse().inverse();

    return pimpl->M;
  }

  const MatrixXd& OrientationFeature::computeProjectedMassInverse(const Feature& featureDes) const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[OrientationFeature::computeProjectedMassInverse] feature does not depend on configuration");

    const MatrixXd& J = computeJacobian(featureDes);
    const MatrixXd& Minv = pimpl->controlFrame->getModel().getInertiaMatrixInverse();
    pimpl->Minv = J * Minv * J.transpose();

    return pimpl->Minv;
  }

  const MatrixXd& OrientationFeature::computeProjectedMassInverse() const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[OrientationFeature::computeProjectedMassInverse] feature does not depend on configuration");

    const MatrixXd& J = computeJacobian();
    const MatrixXd& Minv = pimpl->controlFrame->getModel().getInertiaMatrixInverse();
    pimpl->Minv = J * Minv * J.transpose();

    return pimpl->Minv;
  }
  TaskState OrientationFeature::getState() const
  {
    TaskState state;
#ifndef OCRA_USES_KDL
      state.setPosition(pimpl->controlFrame->getPosition());
      state.setVelocity(pimpl->controlFrame->getVelocity());
      state.setAcceleration(pimpl->controlFrame->getAcceleration());
      state.setWrench(pimpl->controlFrame->getWrench());
#else
      state.setPositionKDL(pimpl->controlFrame->getPositionKDL());
      state.setVelocityKDL(pimpl->controlFrame->getVelocityKDL());
      state.setAccelerationKDL(pimpl->controlFrame->getAccelerationKDL());
      state.setWrenchKDL(pimpl->controlFrame->getWrenchKDL());
#endif
    return state;
  }

  void OrientationFeature::setState(const TaskState& newState)
  {
#ifndef OCRA_USES_KDL
      try {
          TargetFrame::Ptr targetFrame = std::dynamic_pointer_cast<TargetFrame>(pimpl->controlFrame);
          if(newState.hasPosition()) {
                targetFrame->setPosition(newState.getPosition());
            }
            if(newState.hasVelocity()) {
                targetFrame->setVelocity(newState.getVelocity());
            }
            if(newState.hasAcceleration()) {
                targetFrame->setAcceleration(newState.getAcceleration());
            }
            if(newState.hasWrench()) {
                targetFrame->setWrench(newState.getWrench());
            }
      } catch (int errCode) {
          std::cout << "You cannot set the state of this feature because it is not a desired feature. It must be constructed with a TargetFrame." << errCode << std::endl;
      }
#else
      try {
          TargetFrame::Ptr targetFrame = std::dynamic_pointer_cast<TargetFrame>(pimpl->controlFrame);
          if(newState.hasPositionKDL()) {
              targetFrame->setPositionKDL(newState.getPositionKDL());
          }
          if(newState.hasVelocityKDL()) {
              targetFrame->setVelocityKDL(newState.getVelocityKDL());
          }
          if(newState.hasAccelerationKDL()) {
              targetFrame->setAccelerationKDL(newState.getAccelerationKDL());
          }
          if(newState.hasWrenchKDL()) {
              targetFrame->setWrenchKDL(newState.getWrenchKDL());
          }
      } catch (int errCode) {
          std::cout << "You cannot set the state of this feature because it is not a desired feature. It must be constructed with a TargetFrame." << errCode << std::endl;
      }
#endif
  }

  // --- DISPLACEMENT -------------------------------------------

  struct DisplacementFeature::Pimpl
  {
    ControlFrame::Ptr controlFrame;
    ECartesianDof axes;
    int dim;
    MatrixXd u;
    MatrixXd spaceTransform;
    VectorXd error;
    VectorXd errorDot;
    VectorXd effort;
    VectorXd acceleration;
    MatrixXd jacobian;
    MatrixXd M;
    MatrixXd Minv;

    Pimpl(ControlFrame::Ptr cf, ECartesianDof a)
      : controlFrame(cf)
      , dim(3 + utils::computeDimensionFor(axes, NONE))
      , axes(a)
    {
      u.resize(3, dim-3);
      int k = 0;
      if(axes & X)
        u.col(k++) = Vector3d(1., 0., 0.);
      if(axes & Y)
        u.col(k++) = Vector3d(0., 1., 0.);
      if(axes & Z)
        u.col(k++) = Vector3d(0., 0., 1.);

      error = VectorXd::Zero(dim);
      errorDot = VectorXd::Zero(dim);
      effort = VectorXd::Zero(dim);
      acceleration = VectorXd::Zero(dim);
      jacobian = MatrixXd::Zero(dim, cf->getJacobian().cols());
      spaceTransform = MatrixXd(dim, 6);
    }
  };

  DisplacementFeature::DisplacementFeature(const std::string& name, ControlFrame::Ptr frame, ECartesianDof axes)
    : Feature(name)
    , pimpl(new Pimpl(frame, axes))
  {
  }

  const MatrixXd& DisplacementFeature::getSpaceTransform() const
  {
      pimpl->spaceTransform.topRows(3).setIdentity();
      #ifdef OCRA_USES_KDL
          pimpl->spaceTransform.bottomRows(3) = pimpl->u.transpose() * util::getKDLFrameAdjoint(pimpl->controlFrame->getPositionKDL());
      #else
          const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
          pimpl->spaceTransform.bottomRows(3) = pimpl->u.transpose() * R.adjoint();
      #endif


    return pimpl->spaceTransform;
  }

  int DisplacementFeature::getDimension() const
  {
    return pimpl->dim;
  }

  const VectorXd& DisplacementFeature::computeEffort(const Feature& featureDes) const
  {
      const DisplacementFeature& sdes = dynamic_cast<const DisplacementFeature&>(featureDes);

      #ifdef OCRA_USES_KDL
        pimpl->effort = pimpl->u.transpose() * util::KDLWrenchToEigenVectorXd(pimpl->controlFrame->getWrenchKDL() - sdes.pimpl->controlFrame->getWrenchKDL());
      #else

         pimpl->effort = pimpl->u.transpose() * (pimpl->controlFrame->getWrench() - sdes.pimpl->controlFrame->getWrench());

        //   // Twist error in the mobile frame
        //   const Eigen::Displacementd Herror = sdes.pimpl->controlFrame->getPosition().inverse() * pimpl->controlFrame->getPosition();
        //   const Eigen::Wrenchd Werror = pimpl->controlFrame->getWrench() - Herror.adjointTr(sdes.pimpl->controlFrame->getWrench());
          //
        //   // project the translational part on the controlled axes
        //   const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
        //   const MatrixXd u_in_mobileFrame = R.inverse().adjoint() * pimpl->u;
        //   pimpl->effort.tail(pimpl->dim - 3) = u_in_mobileFrame.transpose() * Werror.getForce();
          //
        //   pimpl->effort.head(3) = Werror.getTorque();

      #endif


    return pimpl->effort;
  }

  const VectorXd& DisplacementFeature::computeEffort() const
  {
      #ifdef OCRA_USES_KDL
          pimpl->effort = pimpl->u.transpose() * util::KDLWrenchToEigenVectorXd(pimpl->controlFrame->getWrenchKDL());
      #else
         pimpl->effort = pimpl->u.transpose() * pimpl->controlFrame->getWrench();
        //   // Twist error in the mobile frame
        //   const Eigen::Wrenchd Werror = pimpl->controlFrame->getWrench();
          //
        //   // project the translational part on the controlled axes
        //   const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
        //   const MatrixXd u_in_mobileFrame = R.inverse().adjoint() * pimpl->u;
        //   pimpl->effort.tail(pimpl->dim - 3) = u_in_mobileFrame.transpose() * Werror.getForce();
        //   pimpl->effort.head(3) = Werror.getTorque();
      #endif


    return pimpl->effort;
  }

  const VectorXd& DisplacementFeature::computeAcceleration(const Feature& featureDes) const
  {
    const DisplacementFeature& sdes = dynamic_cast<const DisplacementFeature&>(featureDes);

    #ifdef OCRA_USES_KDL
        pimpl->acceleration = pimpl->u.transpose() * util::KDLTwistToEigenVectorXd(pimpl->controlFrame->getAccelerationKDL() - sdes.pimpl->controlFrame->getAccelerationKDL());
    #else
        pimpl->acceleration = pimpl->u.transpose() * (pimpl->controlFrame->getAcceleration() - sdes.pimpl->controlFrame->getAcceleration());
    #endif


    return pimpl->acceleration;
  }

  const VectorXd& DisplacementFeature::computeAcceleration() const
  {
      #ifdef OCRA_USES_KDL
          pimpl->acceleration = pimpl->u.transpose() * util::KDLTwistToEigenVectorXd(pimpl->controlFrame->getAccelerationKDL());
      #else
          pimpl->acceleration =pimpl->u.transpose() * pimpl->controlFrame->getAcceleration();
      #endif

    return pimpl->acceleration;
  }

const VectorXd& DisplacementFeature::computeError(const Feature& featureDes) const
{
    const DisplacementFeature& sdes = dynamic_cast<const DisplacementFeature&>(featureDes);

    #ifdef OCRA_USES_KDL
        const KDL::Frame tmp = sdes.pimpl->controlFrame->getPositionKDL().Inverse()*pimpl->controlFrame->getPositionKDL();
        Eigen::Vector3d tmpLog;
        util::quaternionLogFromKDLFrame(tmp, tmpLog);
        pimpl->error.head(3) = util::rotationMatrixFromKDLFrame(sdes.pimpl->controlFrame->getPositionKDL())*tmpLog;
        pimpl->error.tail(pimpl->dim - 3) = pimpl->u.transpose() * util::KDLVectorToEigenVector3d(pimpl->controlFrame->getPositionKDL().p - sdes.pimpl->controlFrame->getPositionKDL().p);
    #else
        const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
        const Eigen::Displacementd::Rotation3D& Rdes = sdes.pimpl->controlFrame->getPosition().getRotation();
        pimpl->error.head(3) = Rdes.adjoint()*((Rdes.inverse() * R).log());
        pimpl->error.tail(pimpl->dim - 3) = pimpl->u.transpose() * (pimpl->controlFrame->getPosition().getTranslation() - sdes.pimpl->controlFrame->getPosition().getTranslation());
    #endif

    return pimpl->error;
}

  const VectorXd& DisplacementFeature::computeError() const
  {
      #ifdef OCRA_USES_KDL
          Eigen::Vector3d tmpLog;
          util::quaternionLogFromKDLFrame(pimpl->controlFrame->getPositionKDL(), tmpLog);
           pimpl->error.head(3) = tmpLog;
          pimpl->error.tail(pimpl->dim - 3) = pimpl->u.transpose() * util::KDLVectorToEigenVector3d(pimpl->controlFrame->getPositionKDL().p);
      #else
        // // Displacement error in the mobile frame
        // //    const Eigen::Displacementd Herror = pimpl->controlFrame->getPosition().inverse();
        //
        // // Project the opposite translational part on the controlled axes
        // const Eigen::Displacementd::Rotation3D& R = pimpl->controlFrame->getPosition().getRotation();
        // const MatrixXd u_in_mobileFrame = R.inverse().adjoint() * pimpl->u;
        // //    pimpl->error.tail(pimpl->dim - 3) = - u_in_mobileFrame.transpose() * Herror.getTranslation();
        pimpl->error.tail(pimpl->dim - 3) = pimpl->u.transpose() *  pimpl->controlFrame->getPosition().getTranslation();
        pimpl->error.head(3) = pimpl->controlFrame->getPosition().getRotation().log();

      #endif
    return pimpl->error;
  }

  const VectorXd& DisplacementFeature::computeErrorDot(const Feature& featureDes) const
  {
    const DisplacementFeature& sdes = dynamic_cast<const DisplacementFeature&>(featureDes);

    #ifdef OCRA_USES_KDL
        pimpl->errorDot = util::KDLTwistToEigenVectorXd(pimpl->controlFrame->getVelocityKDL() - sdes.pimpl->controlFrame->getVelocityKDL());

        pimpl->errorDot.tail(pimpl->dim - 3) = (pimpl->u.transpose() * pimpl->errorDot.tail(pimpl->dim - 3)).eval();

    #else
        pimpl->errorDot.tail(pimpl->dim - 3) = pimpl->u.transpose() * (pimpl->controlFrame->getVelocity().getLinearVelocity() - sdes.pimpl->controlFrame->getVelocity().getLinearVelocity());

        pimpl->errorDot.head(3) = pimpl->controlFrame->getVelocity().getAngularVelocity() - sdes.pimpl->controlFrame->getVelocity().getAngularVelocity();
    #endif


    return pimpl->errorDot;
  }

  const VectorXd& DisplacementFeature::computeErrorDot() const
  {
      #ifdef OCRA_USES_KDL
          pimpl->errorDot = util::KDLTwistToEigenVectorXd(pimpl->controlFrame->getVelocityKDL());

          pimpl->errorDot.tail(pimpl->dim - 3) = (pimpl->u.transpose() * pimpl->errorDot.tail(pimpl->dim - 3)).eval();

      #else
          pimpl->errorDot.tail(pimpl->dim - 3) = pimpl->u.transpose() * pimpl->controlFrame->getVelocity().getLinearVelocity();

          pimpl->errorDot.head(3) = pimpl->controlFrame->getVelocity().getAngularVelocity();
      #endif

    return pimpl->errorDot;
  }

  const MatrixXd& DisplacementFeature::computeJacobian(const Feature& featureDes) const
  {
    const DisplacementFeature& sdes = dynamic_cast<const DisplacementFeature&>(featureDes);

    const MatrixXd J = pimpl->controlFrame->getJacobian() - sdes.pimpl->controlFrame->getJacobian();

    return pimpl->jacobian;
  }

  const MatrixXd& DisplacementFeature::computeJacobian() const
  {
    pimpl->jacobian = pimpl->controlFrame->getJacobian();
    return pimpl->jacobian;
  }

  const MatrixXd& DisplacementFeature::computeProjectedMass(const Feature& featureDes) const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[DisplacementFeature::computeProjectedMass] feature does not depend on configuration");

    pimpl->M = computeProjectedMassInverse(featureDes).inverse();

    return pimpl->M;
  }

  const MatrixXd& DisplacementFeature::computeProjectedMass() const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[DisplacementFeature::computeProjectedMass] feature does not depend on configuration");

    pimpl->M = computeProjectedMassInverse().inverse();

    return pimpl->M;
  }

  const MatrixXd& DisplacementFeature::computeProjectedMassInverse(const Feature& featureDes) const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[DisplacementFeature::computeProjectedMassInverse] feature does not depend on configuration");

    const MatrixXd& J = computeJacobian(featureDes);
    const MatrixXd& Minv = pimpl->controlFrame->getModel().getInertiaMatrixInverse();
    pimpl->Minv = J * Minv * J.transpose();

    return pimpl->Minv;
  }

  const MatrixXd& DisplacementFeature::computeProjectedMassInverse() const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[DisplacementFeature::computeProjectedMassInverse] feature does not depend on configuration");

    const MatrixXd& J = computeJacobian();
    const MatrixXd& Minv = pimpl->controlFrame->getModel().getInertiaMatrix().inverse();
    pimpl->Minv = J * Minv * J.transpose();

    return pimpl->Minv;
  }
  TaskState DisplacementFeature::getState() const
  {
      TaskState state;
      #ifdef OCRA_USES_KDL
        state.setPositionKDL(pimpl->controlFrame->getPositionKDL());
        state.setVelocityKDL(pimpl->controlFrame->getVelocityKDL());
        state.setAccelerationKDL(pimpl->controlFrame->getAccelerationKDL());
        state.setWrenchKDL(pimpl->controlFrame->getWrenchKDL());
      #else
          state.setPosition(pimpl->controlFrame->getPosition());
          state.setVelocity(pimpl->controlFrame->getVelocity());
          state.setAcceleration(pimpl->controlFrame->getAcceleration());
          state.setWrench(pimpl->controlFrame->getWrench());
      #endif

      return state;
  }

  void DisplacementFeature::setState(const TaskState& newState)
  {
      #ifdef OCRA_USES_KDL
      try {
          TargetFrame::Ptr targetFrame = std::dynamic_pointer_cast<TargetFrame>(pimpl->controlFrame);
          if(newState.hasPositionKDL()) {
                targetFrame->setPositionKDL(newState.getPositionKDL());
            }
            if(newState.hasVelocityKDL()) {
                targetFrame->setVelocityKDL(newState.getVelocityKDL());
            }
            if(newState.hasAccelerationKDL()) {
                targetFrame->setAccelerationKDL(newState.getAccelerationKDL());
            }
            if(newState.hasWrenchKDL()) {
                targetFrame->setWrenchKDL(newState.getWrenchKDL());
            }
      } catch (int errCode) {
          std::cout << "You cannot set the state of this feature because it is not a desired feature. It must be constructed with a TargetFrame." << errCode << std::endl;
      }
      #else
          try {
              TargetFrame::Ptr targetFrame = std::dynamic_pointer_cast<TargetFrame>(pimpl->controlFrame);
              if(newState.hasPosition()) {
                    targetFrame->setPosition(newState.getPosition());
                }
                if(newState.hasVelocity()) {
                    targetFrame->setVelocity(newState.getVelocity());
                }
                if(newState.hasAcceleration()) {
                    targetFrame->setAcceleration(newState.getAcceleration());
                }
                if(newState.hasWrench()) {
                    targetFrame->setWrench(newState.getWrench());
                }
          } catch (int errCode) {
              std::cout << "You cannot set the state of this feature because it is not a desired feature. It must be constructed with a TargetFrame." << errCode << std::endl;
          }
      #endif


  }

  // --- CONTACT CONSTRAINT FEATURES ----------------------------

  struct ContactConstraintFeature::Pimpl
  {
    ControlFrame::Ptr controlFrame;
    MatrixXd spaceTransform;
    VectorXd error;
    VectorXd errorDot;
    VectorXd effort;
    VectorXd acceleration;
    MatrixXd jacobian;
    MatrixXd M;
    MatrixXd Minv;

    Pimpl(ControlFrame::Ptr cf)
      : controlFrame(cf)
      , spaceTransform(MatrixXd::Identity(6, 6))
      , error(VectorXd::Zero(6))
      , errorDot(VectorXd::Zero(6))
      , effort(VectorXd::Zero(6))
      , acceleration(VectorXd::Zero(6))
      , jacobian(MatrixXd::Zero(6, cf->getJacobian().cols()))
    {
    }
  };

  ContactConstraintFeature::ContactConstraintFeature(const std::string& name, ControlFrame::Ptr frame)
    : Feature(name)
    , pimpl(new Pimpl(frame))
  {
  }

  const MatrixXd& ContactConstraintFeature::getSpaceTransform() const
  {
    return pimpl->spaceTransform;
  }

  int ContactConstraintFeature::getDimension() const
  {
    return 6;
  }

  const VectorXd& ContactConstraintFeature::computeEffort(const Feature& featureDes) const
  {
    throw std::runtime_error("[ContactConstraintFeature::computeEffort(const Feature&)] Desired feature are irrelevant in ContactConstraintFeature");
  }

  const VectorXd& ContactConstraintFeature::computeEffort() const
  {
    const Eigen::Wrenchd Werror = pimpl->controlFrame->getWrench();
    pimpl->effort.tail(3) = Werror.getForce();
    pimpl->effort.head(3) = Werror.getTorque();
    return pimpl->effort;
  }

  const VectorXd& ContactConstraintFeature::computeAcceleration(const Feature& featureDes) const
  {
    throw std::runtime_error("[ContactConstraintFeature::computeAcceleration(const Feature&)] Desired feature are irrelevant in ContactConstraintFeature");
  }

  const VectorXd& ContactConstraintFeature::computeAcceleration() const
  {
    pimpl->acceleration.head(3) = pimpl->controlFrame->getAcceleration().getAngularVelocity();
    pimpl->acceleration.tail(3) = pimpl->controlFrame->getAcceleration().getLinearVelocity();
    return pimpl->acceleration;
  }

  const VectorXd& ContactConstraintFeature::computeError(const Feature& featureDes) const
  {
    throw std::runtime_error("[ContactConstraintFeature::computeError(const Feature&)] Desired feature are irrelevant in ContactConstraintFeature");
  }

  const VectorXd& ContactConstraintFeature::computeError() const
  {
    const Eigen::Displacementd H = pimpl->controlFrame->getPosition();
    const Eigen::Displacementd::Rotation3D& R = H.getRotation();
//    pimpl->error.tail(3) = R.adjoint().transpose() * H.getTranslation();
    pimpl->error.tail(3) = H.getTranslation();
    pimpl->error.head(3) = R.log();
    return pimpl->error;
  }

  const VectorXd& ContactConstraintFeature::computeErrorDot(const Feature& featureDes) const
  {
    throw std::runtime_error("[ContactConstraintFeature::computeErrorDot(const Feature&)] Desired feature are irrelevant in ContactConstraintFeature");
  }

  const VectorXd& ContactConstraintFeature::computeErrorDot() const
  {
    const Eigen::Twistd Terror = pimpl->controlFrame->getVelocity();
    pimpl->errorDot.tail(3) = Terror.getLinearVelocity();
    pimpl->errorDot.head(3) = Terror.getAngularVelocity();
    return pimpl->errorDot;
  }

  const MatrixXd& ContactConstraintFeature::computeJacobian(const Feature& featureDes) const
  {
    throw std::runtime_error("[ContactConstraintFeature::computeJacobian(const Feature&)] Desired feature are irrelevant in ContactConstraintFeature");
  }

  const MatrixXd& ContactConstraintFeature::computeJacobian() const
  {
    pimpl->jacobian = pimpl->controlFrame->getJacobian();
    return pimpl->jacobian;
  }

  const MatrixXd& ContactConstraintFeature::computeProjectedMass(const Feature& featureDes) const
  {
    throw std::runtime_error("[ContactConstraintFeature::computeProjectedMass(const Feature&)] Desired feature are irrelevant in ContactConstraintFeature");
  }

  const MatrixXd& ContactConstraintFeature::computeProjectedMass() const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[DisplacementFeature::computeProjectedMass] feature does not depend on configuration");

    pimpl->M = computeProjectedMassInverse().inverse();

    return pimpl->M;
  }

  const MatrixXd& ContactConstraintFeature::computeProjectedMassInverse(const Feature& featureDes) const
  {
    throw std::runtime_error("[ContactConstraintFeature::computeProjectedMassInverse(const Feature&)] Desired feature are irrelevant in ContactConstraintFeature");
  }

  const MatrixXd& ContactConstraintFeature::computeProjectedMassInverse() const
  {
    if(!pimpl->controlFrame->dependsOnModelConfiguration())
      throw std::runtime_error("[DisplacementFeature::computeProjectedMassInverse] feature does not depend on configuration");

    const MatrixXd& J = computeJacobian();
    const MatrixXd& Minv = pimpl->controlFrame->getModel().getInertiaMatrix().inverse();
    pimpl->Minv = J * Minv * J.transpose();

    return pimpl->Minv;
  }
  TaskState ContactConstraintFeature::getState() const
  {
      TaskState state;
      state.setPosition(pimpl->controlFrame->getPosition());
      state.setVelocity(pimpl->controlFrame->getVelocity());
      state.setAcceleration(pimpl->controlFrame->getAcceleration());
      state.setWrench(pimpl->controlFrame->getWrench());

      return state;
  }

  void ContactConstraintFeature::setState(const TaskState& newState)
  {
      try {
          TargetFrame::Ptr targetFrame = std::dynamic_pointer_cast<TargetFrame>(pimpl->controlFrame);
          if(newState.hasPosition()) {
                targetFrame->setPosition(newState.getPosition());
            }
            if(newState.hasVelocity()) {
                targetFrame->setVelocity(newState.getVelocity());
            }
            if(newState.hasAcceleration()) {
                targetFrame->setAcceleration(newState.getAcceleration());
            }
            if(newState.hasWrench()) {
                targetFrame->setWrench(newState.getWrench());
            }
      } catch (int errCode) {
          std::cout << "You cannot set the state of this feature because it is not a desired feature. It must be constructed with a TargetFrame." << errCode << std::endl;
      }
  }




  // --- ARTICULAR ----------------------------------------------

  struct FullStateFeature::Pimpl
  {
    FullState::Ptr state;
    VectorXd error;
    VectorXd errorDot;
    VectorXd effort;
    VectorXd acceleration;
    MatrixXd J;
    MatrixXd M;
    MatrixXd Minv;
    MatrixXd spaceTransform;

    Pimpl(FullState::Ptr fs)
      : state(fs)
    {
      spaceTransform = MatrixXd::Identity(state->getSize(), state->getSize());
    }
  };

  FullStateFeature::FullStateFeature(const std::string& name, FullState::Ptr state)
    : Feature(name)
    , pimpl( new Pimpl(state) )
  {
  }

  const MatrixXd& FullStateFeature::getSpaceTransform() const
  {
    return pimpl->spaceTransform;
  }

  int FullStateFeature::getDimension() const
  {
    return pimpl->state->getSize();
  }

  const VectorXd& FullStateFeature::computeEffort(const Feature& featureDes) const
  {
    const FullStateFeature& sdes = dynamic_cast<const FullStateFeature&>(featureDes);
    pimpl->effort = pimpl->state->tau() - sdes.pimpl->state->tau();
    return pimpl->effort;
  }

  const VectorXd& FullStateFeature::computeEffort() const
  {
    pimpl->effort = pimpl->state->tau();
    return pimpl->effort;
  }

  const VectorXd& FullStateFeature::computeAcceleration(const Feature& featureDes) const
  {
    const FullStateFeature& sdes = dynamic_cast<const FullStateFeature&>(featureDes);
    pimpl->acceleration = pimpl->state->qddot() - sdes.pimpl->state->qddot();
    return pimpl->acceleration;
  }

  const VectorXd& FullStateFeature::computeAcceleration() const
  {
    pimpl->acceleration = pimpl->state->qddot();
    return pimpl->acceleration;
  }

  const VectorXd& FullStateFeature::computeError(const Feature& featureDes) const
  {
    const FullStateFeature& sdes = dynamic_cast<const FullStateFeature&>(featureDes);
    pimpl->error = pimpl->state->q() - sdes.pimpl->state->q();
    return pimpl->error;
  }

  const VectorXd& FullStateFeature::computeError() const
  {
    pimpl->error = pimpl->state->q();
    return pimpl->error;
  }

  const VectorXd& FullStateFeature::computeErrorDot(const Feature& featureDes) const
  {
    const FullStateFeature& sdes = dynamic_cast<const FullStateFeature&>(featureDes);
    pimpl->errorDot = pimpl->state->qdot() - sdes.pimpl->state->qdot();
    return pimpl->errorDot;
  }

  const VectorXd& FullStateFeature::computeErrorDot() const
  {
    pimpl->errorDot = pimpl->state->qdot();
    return pimpl->errorDot;
  }

  const MatrixXd& FullStateFeature::computeJacobian(const Feature& featureDes) const
  {
    const FullStateFeature& sdes = dynamic_cast<const FullStateFeature&>(featureDes);
    pimpl->J = pimpl->state->getJacobian();
    return pimpl->J;
  }

  const MatrixXd& FullStateFeature::computeJacobian() const
  {
    pimpl->J = pimpl->state->getJacobian();
    return pimpl->J;
  }

  const MatrixXd& FullStateFeature::computeProjectedMass(const Feature& featureDes) const
  {
    const FullStateFeature& sdes = dynamic_cast<const FullStateFeature&>(featureDes);
    pimpl->M = pimpl->state->getInertiaMatrix();
    return pimpl->M;
  }

  const MatrixXd& FullStateFeature::computeProjectedMass() const
  {
    pimpl->M = pimpl->state->getInertiaMatrix();
    return pimpl->M;
  }

  const MatrixXd& FullStateFeature::computeProjectedMassInverse(const Feature& featureDes) const
  {
    const FullStateFeature& sdes = dynamic_cast<const FullStateFeature&>(featureDes);
    pimpl->M = pimpl->state->getInertiaMatrixInverse();
    return pimpl->M;
  }

  const MatrixXd& FullStateFeature::computeProjectedMassInverse() const
  {
    pimpl->M = pimpl->state->getInertiaMatrixInverse();
    return pimpl->M;
  }

  TaskState FullStateFeature::getState() const
  {
      TaskState state;
      state.setQ(pimpl->state->q());
      state.setQd(pimpl->state->qdot());
      state.setQdd(pimpl->state->qddot());

      return state;
  }

  void FullStateFeature::setState(const TaskState& newState)
  {
      try {
          FullTargetState::Ptr targetState = std::dynamic_pointer_cast<FullTargetState>(pimpl->state);
          if(newState.hasQ()) {
                targetState->set_q(newState.getQ());
            }
            if(newState.hasQd()) {
                targetState->set_qdot(newState.getQd());
            }
            if(newState.hasQdd()) {
                targetState->set_qddot(newState.getQdd());
            }
      } catch (int errCode) {
          std::cout << "You cannot set the state of this feature because it is not a desired feature. It must be constructed with a FullTargetState." << errCode << std::endl;
      }
  }

  // --- PARTIAL - ARTICULAR ------------------------------------

  struct PartialStateFeature::Pimpl
  {
      PartialState::Ptr state;
      VectorXd error;
      VectorXd errorDot;
      VectorXd effort;
      VectorXd acceleration;
      MatrixXd J;
      MatrixXd M;
      MatrixXd Minv;
      MatrixXd spaceTransform;

      Pimpl(PartialState::Ptr ps)
          : state(ps)
      {
          spaceTransform = MatrixXd::Identity(state->getSize(), state->getSize());
      }
    };

  PartialStateFeature::PartialStateFeature(const std::string& name, PartialState::Ptr state)
      : Feature(name)
      , pimpl( new Pimpl(state) )
  {
  }

  const MatrixXd& PartialStateFeature::getSpaceTransform() const
  {
      return pimpl->spaceTransform;
  }

  int PartialStateFeature::getDimension() const
  {
      return pimpl->state->getSize();
  }

  const VectorXd& PartialStateFeature::computeEffort(const Feature& featureDes) const
  {
      const PartialStateFeature& sdes = dynamic_cast<const PartialStateFeature&>(featureDes);
      pimpl->effort = pimpl->state->tau() - sdes.pimpl->state->tau();
      return pimpl->effort;
  }

  const VectorXd& PartialStateFeature::computeEffort() const
  {
      pimpl->effort = pimpl->state->tau();
      return pimpl->effort;
  }

  const VectorXd& PartialStateFeature::computeAcceleration(const Feature& featureDes) const
  {
      const PartialStateFeature& sdes = dynamic_cast<const PartialStateFeature&>(featureDes);
      pimpl->acceleration = pimpl->state->qddot() - sdes.pimpl->state->qddot();
      return pimpl->acceleration;
  }

  const VectorXd& PartialStateFeature::computeAcceleration() const
  {
      pimpl->acceleration = pimpl->state->qddot();
      return pimpl->acceleration;
  }

  const VectorXd& PartialStateFeature::computeError(const Feature& featureDes) const
  {
      const PartialStateFeature& sdes = dynamic_cast<const PartialStateFeature&>(featureDes);
      pimpl->error = pimpl->state->q() - sdes.pimpl->state->q();
      return pimpl->error;
  }

  const VectorXd& PartialStateFeature::computeError() const
  {
      pimpl->error = pimpl->state->q();
      return pimpl->error;
  }

  const VectorXd& PartialStateFeature::computeErrorDot(const Feature& featureDes) const
  {
      const PartialStateFeature& sdes = dynamic_cast<const PartialStateFeature&>(featureDes);
      pimpl->errorDot = pimpl->state->qdot() - sdes.pimpl->state->qdot();
      return pimpl->errorDot;
  }

  const VectorXd& PartialStateFeature::computeErrorDot() const
  {
      pimpl->errorDot = pimpl->state->qdot();
      return pimpl->errorDot;
  }

  const MatrixXd& PartialStateFeature::computeJacobian(const Feature& featureDes) const
  {
      const PartialStateFeature& sdes = dynamic_cast<const PartialStateFeature&>(featureDes);
      pimpl->J = pimpl->state->getJacobian();
      return pimpl->J;
  }

  const MatrixXd& PartialStateFeature::computeJacobian() const
  {
      pimpl->J = pimpl->state->getJacobian();
      return pimpl->J;
  }

  const MatrixXd& PartialStateFeature::computeProjectedMass(const Feature& featureDes) const
  {
      const PartialStateFeature& sdes = dynamic_cast<const PartialStateFeature&>(featureDes);
      pimpl->M = pimpl->state->getInertiaMatrix();
      return pimpl->M;
  }

  const MatrixXd& PartialStateFeature::computeProjectedMass() const
  {
      pimpl->M = pimpl->state->getInertiaMatrix();
      return pimpl->M;
  }

  const MatrixXd& PartialStateFeature::computeProjectedMassInverse(const Feature& featureDes) const
  {
      const PartialStateFeature& sdes = dynamic_cast<const PartialStateFeature&>(featureDes);
      pimpl->M = pimpl->state->getInertiaMatrixInverse();
      return pimpl->M;
  }

  const MatrixXd& PartialStateFeature::computeProjectedMassInverse() const
  {
      pimpl->M = pimpl->state->getInertiaMatrixInverse();
      return pimpl->M;
  }
  TaskState PartialStateFeature::getState() const
  {
      TaskState state;
      state.setQ(pimpl->state->q());
      state.setQd(pimpl->state->qdot());
      state.setQdd(pimpl->state->qddot());

      return state;
  }

  void PartialStateFeature::setState(const TaskState& newState)
  {
      try {
          PartialTargetState::Ptr targetState = std::dynamic_pointer_cast<PartialTargetState>(pimpl->state);
          if(newState.hasQ()) {
                targetState->set_q(newState.getQ());
            }
            if(newState.hasQd()) {
                targetState->set_qdot(newState.getQd());
            }
            if(newState.hasQdd()) {
                targetState->set_qddot(newState.getQdd());
            }
      } catch (int errCode) {
          std::cout << "You cannot set the state of this feature because it is not a desired feature. It must be constructed with a PartialTargetState." << errCode << std::endl;
      }
  }

}

// cmake:sourcegroup=Api
