#include "ocra/control/ControlFrame.h"

#include "ocra/control/Model.h"
#include <iostream>

namespace ocra
{
  // --- ABSTRACT -----------------------------------------------

  ControlFrame::ControlFrame(const std::string& name)
    : NamedInstance(name)
  {
  }

  ControlFrame::~ControlFrame()
  {
  }


  // --- EXTERNAL INPUT -----------------------------------------

  struct TargetFrame::Pimpl
  {
    const Model& model;
    Eigen::Displacementd H;
    Eigen::Twistd T;
    Eigen::Twistd gamma;
    Eigen::Wrenchd W;
    Jacobian6d J;
      //KDL Migration
      KDL::Frame H_KDL;
      KDL::Twist T_KDL;
      KDL::Twist gamma_KDL;

    Pimpl(const Model& m)
      : model(m)
      , H(Displacementd::Identity())
      , T(Twistd::Zero())
      , gamma(Twistd::Zero())
      , W(Wrenchd::Zero())
      , J(MatrixXd::Zero(6, m.nbDofs()))
      , H_KDL(KDL::Frame::Identity())
      , T_KDL(KDL::Twist::Zero())
      , gamma_KDL(KDL::Twist::Zero())
    {}
  };

  TargetFrame::TargetFrame(const std::string& name, const Model& model)
    : ControlFrame(name)
    , pimpl(new Pimpl(model))
  {
  }

  Eigen::Displacementd TargetFrame::getPosition() const
  {
    return pimpl->H;
  }

  Eigen::Twistd TargetFrame::getVelocity() const
  {
    return pimpl->T;
  }

  Eigen::Twistd TargetFrame::getAcceleration() const
  {
    return pimpl->gamma;
  }

  Eigen::Wrenchd TargetFrame::getWrench() const
  {
    return pimpl->W;
  }

  Jacobian6d TargetFrame::getJacobian() const
  {
    return pimpl->J;
  }

  bool TargetFrame::dependsOnModelConfiguration() const
  {
    return false;
  }

  const Model& TargetFrame::getModel() const
  {
    return pimpl->model;
  }

  void TargetFrame::setPosition(const Eigen::Displacementd& H)
  {
    pimpl->H = H;
  }

  void TargetFrame::setVelocity(const Eigen::Twistd& T)
  {
    pimpl->T = T;
  }

  void TargetFrame::setAcceleration(const Eigen::Twistd& gamma)
  {
    pimpl->gamma = gamma;
  }

  void TargetFrame::setWrench(const Eigen::Wrenchd& W)
  {
    pimpl->W = W;
  }
    
    // KDL MIGRATION
    KDL::Frame TargetFrame::getPositionKDL() const
    {
        return pimpl->H_KDL;
    }
    
    KDL::Twist TargetFrame::getVelocityKDL() const
    {
        return pimpl->T_KDL;
    }
    
    KDL::Twist TargetFrame::getAccelerationKDL() const
    {
        return pimpl->gamma_KDL;
    }
    
    void TargetFrame::setPositionKDL(const KDL::Frame &H)
    {
        pimpl->H_KDL = H;
    }
    
    void TargetFrame::setVelocityKDL(const KDL::Twist &T)
    {
        pimpl->T_KDL = T;
    }
    
    void TargetFrame::setAccelerationKDL(const KDL::Twist &gamma)
    {
        pimpl->gamma_KDL = gamma;
    }


  // --- ATTACHED TO A SEGMENT ----------------------------------

  struct SegmentFrame::Pimpl
  {
    const Model& model;
    int index;
    Eigen::Displacementd H_localFrame;
    MatrixXd Adj_H_segment_in_controlledFrame; // TODO [minor]: fixed size matrix
      // KDL migration
      KDL::Frame H_localFrameKDL;

    Pimpl(const Model& m, const std::string& segname)
      : model(m)
      , index(model.getSegmentIndex(segname))
      , H_localFrame(Displacementd::Identity())
      , Adj_H_segment_in_controlledFrame(MatrixXd::Identity(6, 6))
      , H_localFrameKDL(KDL::Frame::Identity())
    {}

    Pimpl(const Model& m, const std::string& segname, const Eigen::Displacementd& H_local)
      : model(m)
      , index(model.getSegmentIndex(segname))
      , H_localFrame(H_local)
      , Adj_H_segment_in_controlledFrame(H_local.inverse().adjoint())
    {}
      
      // KDL Migration
      Pimpl(const Model& m, const std::string& segname, KDL::Frame& H_local)
      : model(m)
      , index(model.getSegmentIndex(segname))
      , Adj_H_segment_in_controlledFrame(ocra::util::getKDLFrameAdjoint(H_local.Inverse()))
      , H_localFrameKDL(H_local)
      {}


    Pimpl(const Model& m, int segmentId)
      : model(m)
      , index(segmentId)
      , H_localFrame(Displacementd::Identity())
      , Adj_H_segment_in_controlledFrame(MatrixXd::Identity(6, 6))
    {}

    Pimpl(const Model& m, int segmentId, const Eigen::Displacementd& H_local)
      : model(m)
      , index(segmentId)
      , H_localFrame(H_local)
      , Adj_H_segment_in_controlledFrame(H_local.inverse().adjoint())
    {}
      
      // KDL Migration
      //TODO: How do we substitute Adj_H_segment_in_controlledFrame with its KDL version
      Pimpl(const Model& m, int segmentId, const KDL::Frame& H_local)
      : model(m)
      , index(segmentId)
      , H_localFrameKDL(H_local)
      //      , Adj_H_segment_in_controlledFrame(H_local.inverse().adjoint())
      {}

  };

  SegmentFrame::SegmentFrame(const std::string& name, const Model& model, const std::string& segname)
    : ControlFrame(name)
    , pimpl(new Pimpl(model, segname))
  {
  }

  SegmentFrame::SegmentFrame(const std::string& name, const Model& model, const std::string& segname, const Eigen::Displacementd& H_local)
    : ControlFrame(name)
    , pimpl(new Pimpl(model, segname, H_local))
  {
  }

  SegmentFrame::SegmentFrame(const std::string& name, const Model& model, int segmentId)
    : ControlFrame(name)
    , pimpl(new Pimpl(model, segmentId))
  {
  }

  SegmentFrame::SegmentFrame(const std::string& name, const Model& model, int segmentId, const Eigen::Displacementd& H_local)
    : ControlFrame(name)
    , pimpl(new Pimpl(model, segmentId, H_local))
  {
  }

  Eigen::Displacementd SegmentFrame::getPosition() const
  {
    return pimpl->model.getSegmentPosition(pimpl->index) * pimpl->H_localFrame;
  }

  Eigen::Twistd SegmentFrame::getVelocity() const
  {
    return pimpl->Adj_H_segment_in_controlledFrame * pimpl->model.getSegmentVelocity(pimpl->index);
  }

  Eigen::Twistd SegmentFrame::getAcceleration() const
  {
    return Eigen::Twistd::Zero();
  }

  Eigen::Wrenchd SegmentFrame::getWrench() const
  {
    return Eigen::Wrenchd::Zero();
  }

  Jacobian6d SegmentFrame::getJacobian() const
  {
    return pimpl->Adj_H_segment_in_controlledFrame * pimpl->model.getSegmentJacobian(pimpl->index);
  }

  bool SegmentFrame::dependsOnModelConfiguration() const
  {
    return true;
  }

  const Model& SegmentFrame::getModel() const
  {
    return pimpl->model;
  }

  int SegmentFrame::getSegmentIndex() const
  {
    return pimpl->index;
  }


  // --- ATTACHED TO THE COM ------------------------------------

  struct CoMFrame::Pimpl
  {
    const Model& model;

    Pimpl(const Model& m): model(m) {}
  };

  CoMFrame::CoMFrame(const std::string& name, const Model& model)
    : ControlFrame(name)
    , pimpl(new Pimpl(model))
  {
  }

  Eigen::Displacementd CoMFrame::getPosition() const
  {
    Eigen::Vector3d tmp = pimpl->model.getCoMPosition();
    return Eigen::Displacementd(tmp, Quaterniond::Identity());
  }

  Eigen::Twistd CoMFrame::getVelocity() const
  {
    return Eigen::Twistd(Twistd::AngularVelocity(0., 0., 0.), pimpl->model.getCoMVelocity());
  }

  Eigen::Twistd CoMFrame::getAcceleration() const
  {
    return Eigen::Twistd::Zero();
  }

  Eigen::Wrenchd CoMFrame::getWrench() const
  {
    return Eigen::Wrenchd::Zero();
  }

  Jacobian6d CoMFrame::getJacobian() const
  {
    Jacobian6d result(6, pimpl->model.nbDofs());
    result.setZero();
    result.block(3, 0, 3, pimpl->model.nbDofs()) = pimpl->model.getCoMJacobian();
    return result;
  }

  bool CoMFrame::dependsOnModelConfiguration() const
  {
    return true;
  }

  const Model& CoMFrame::getModel() const
  {
    return pimpl->model;
  }
}

// cmake:sourcegroup=Api
