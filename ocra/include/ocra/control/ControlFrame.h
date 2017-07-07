/*!
\file ControlFrame.h
\brief Classes that represent frames.

Copyright (C) 2010 CEA/DRT/LIST/DTSI/SRCI

\author Evrard Paul
\author Escande Adrien
\date 2010/10/03

File history:
*/

#ifndef _OCRA_CONTROL_FRAME_H_
#define _OCRA_CONTROL_FRAME_H_

#include "ocra/optim/NamedInstance.h"
#include <ocra/util/MathTypes.h>
#include <ocra/util/Macros.h>
#include "Eigen/Lgsm"
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <yarp/os/Semaphore.h>

#include "kdl/frames_io.hpp"
#include "kdl/frames.hpp"
#include <ocra/util/ErrorsHelper.h>
#include <ocra/util/KDLUtilities.h>


namespace ocra
{
  class Model;
}

namespace ocra
{
    /**
     * @class ControlFrame
     * @brief Generic representation of a frame. Gives access to its position, velocity, jacobian, etc.
     * @details A frame is always associated to a manikin's model. However, it does not always depend
     on its configuration.
     If it depends on the configuration of the manikin's model, then the method
     ControlFrame::dependsOnModelConfiguration will return true.
     Otherwise, this method returns false, and ControlFrame::getJacobian
     will return a null matrix, whose size be 6 x getModel().nbDofs().
     */
  class ControlFrame
    : public NamedInstance
    , boost::noncopyable
  {
      DEFINE_CLASS_POINTER_TYPEDEFS(ControlFrame)

  protected:
    ControlFrame(const std::string& name);

  public:
    virtual ~ControlFrame() = 0;

      
    /**
     Retrieves the position of this frame w.r.t. the world reference frame.

     @return Frame position.
     */
    virtual Eigen::Displacementd getPosition() const = 0;
      
    /**
     Retrieves the velocity of this frame w.r.t. the world reference frame.

     @return Frame velocity.
     */
    virtual Eigen::Twistd getVelocity() const = 0;
      
    /**
     Returns the acceleration of the frame.

     @return Frame acceleration.
     */
    virtual Eigen::Twistd getAcceleration() const = 0;
      
    /**
     Retrieves the wrench acting on this frame.

     @return Wrench on the frame.
     */
    virtual Eigen::Wrenchd getWrench() const = 0;
      
    /**
     Retrieves the 6xDoF jacobian of the frame, where the angular part is actually
     null. The linear part is retrieved from the model referenced by this class.
     
     @return 6-dimensional frame Jacobian.
     */
    virtual Eigen::Matrix<double,6,Eigen::Dynamic> getJacobian() const = 0;
      
    /**
     Indicates whether the frame depends on the mannequein's configuration.

     @return True when it does, false otherwise.
     */
    virtual bool dependsOnModelConfiguration() const = 0;
      
    /**
     Retrieves a reference to the model passed to this frame.

     @return Reference to model used by this frame.
     */
    virtual const Model& getModel() const = 0;
      
    // ################################# KDL migration ####################################
    /**
    Retrieves the position of this frame w.r.t. the world reference frame.

    @return Frame position.
    */
    virtual KDL::Frame getPositionKDL() const = 0;
      
    /**
    Retrieves the velocity of this frame w.r.t. the world reference frame.

    @return Frame velocity.
    */
    virtual KDL::Twist getVelocityKDL() const = 0;
      
    /**
    Returns the acceleration of the frame.

    @return Frame acceleration.
    */
    virtual KDL::Twist getAccelerationKDL() const = 0;

    /**
    Retrieves the wrench acting on this frame.

    @return Wrench on the frame.
    */
    virtual KDL::Wrench getWrenchKDL() const = 0;
      
  };


  // --- EXTERNAL INPUT -----------------------------------------

    /**
     * @class TargetFrame
     * @brief Represents a 'target', i.e. something that does not depend on a model but must match with another frame.
     * @details A target frame has its position and velocity manually defined via TargetFrame::setPosition and TargetFrame::setVelocity. It can be used to specify control objectives. Note that you always have to specify a displacement and a twist, whatever the control objective. If you are only interested in specifying a position, you can just fill the rotational part with the identity quaternion. Likewise, if only the rotational part is of interest, you can specify any position, since it will be ignored by the controller. The same goes for the twist.
     */

  class TargetFrame
    : public ControlFrame
  {
      DEFINE_CLASS_POINTER_TYPEDEFS(TargetFrame)
  public:
    TargetFrame(const std::string& name, const Model& model);

    Eigen::Displacementd getPosition() const;
    Eigen::Twistd getVelocity() const;
    Eigen::Twistd getAcceleration() const;
    Eigen::Wrenchd getWrench() const;
    Eigen::Matrix<double,6,Eigen::Dynamic> getJacobian() const;
    bool dependsOnModelConfiguration() const;
    const Model& getModel() const;

      
    /**
     Sets the pose of the target frame.

     @param H Pose to be set.
     */
    void setPosition(const Eigen::Displacementd& H);
      
    /**
     Sets the velocity of the target frame.

     @param T Velocity to be set.
     */
    void setVelocity(const Eigen::Twistd& T);
      
    /**
     Sets the acceleration of the target frame.

     @param gamma Acceleration to be set.
     */
    void setAcceleration(const Eigen::Twistd& gamma);
      
    /**
     Sets the wrench at the target frame.

     @param W Wrench to be set.
     */
    void setWrench(const Eigen::Wrenchd& W);
      
      // KDL Changes
      KDL::Frame getPositionKDL() const;
      KDL::Twist getVelocityKDL() const;
      KDL::Twist getAccelerationKDL() const;
      KDL::Wrench getWrenchKDL() const;
      
      /**
       Sets the pose of the target frame.
       
       @param H Pose to be set.
       */
      void setPositionKDL(const KDL::Frame &H);

      /**
       Sets the velocity of the target frame.
       
       @param T Velocity to be set.
       */
      void setVelocityKDL(const KDL::Twist &T);
      
      /**
       Sets the acceleration of the target frame.
       
       @param gamma Acceleration to be set.
       */
      void setAccelerationKDL(const KDL::Twist &gamma);
      
      /**
       Sets the wrench at the target frame.
       
       @param W Wrench to be set.
       */
      void setWrenchKDL(const KDL::Wrench &W);

  private:
    struct Pimpl;
    boost::shared_ptr<Pimpl> pimpl;
  };


  // --- ATTACHED TO A SEGMENT ----------------------------------
    /**
     * @class SegmentFrame
     * @brief  A frame attached to a segment of a model.
     */
  class SegmentFrame
    : public ControlFrame
  {
      DEFINE_CLASS_POINTER_TYPEDEFS(SegmentFrame)
  public:
    SegmentFrame(const std::string& name, const Model& model, const std::string& segname);
    SegmentFrame(const std::string& name, const Model& model, const std::string& segname, const Eigen::Displacementd& H_local);
      // KDL migration of method above
      SegmentFrame(const std::string& name, const Model& model, const std::string& segname, const KDL::Frame& H_local);
    SegmentFrame(const std::string& name, const Model& model, int segmentId);
    SegmentFrame(const std::string& name, const Model& model, int segmentId, const Eigen::Displacementd& H_local);
    // KDL Migration of method above
      SegmentFrame(const std::string& name, const Model& model, int segmentId, const KDL::Frame& H_localKDL);

    Eigen::Displacementd getPosition() const;
    Eigen::Twistd getVelocity() const;
    Eigen::Twistd getAcceleration() const;
    Eigen::Wrenchd getWrench() const;
    Eigen::Matrix<double,6,Eigen::Dynamic> getJacobian() const;
    bool dependsOnModelConfiguration() const;
    const Model& getModel() const;
    int getSegmentIndex() const;
      
      // KDL migration
      KDL::Frame getPositionKDL() const;
      KDL::Twist getVelocityKDL() const;
      KDL::Twist getAccelerationKDL() const;
      KDL::Wrench getWrenchKDL() const;
      
  private:
    struct Pimpl;
    boost::shared_ptr<Pimpl> pimpl;
  };


  // --- ATTACHED TO THE COM ------------------------------------
    /**
     * @class CoMFrame
     * @brief Creates a frame whose center is at the CoM of the model and whose axes are parallel to the axes of the world frame.
     */
  class CoMFrame
    : public ControlFrame
  {
      DEFINE_CLASS_POINTER_TYPEDEFS(CoMFrame)
  public:
    CoMFrame(const std::string& name, const Model& model);
      
    /**
     Returns the position of the CoM as stored in the model object passed to this class. 
     The orientation is set to the identity.

     @return Position of the CoM as a Displacementd object.
     */
    Eigen::Displacementd getPosition() const;
      
    /**
     Retrieves the velocity of the CoM, composed by the linear velocity as stored in the model 
     object passed to this class and null angular velocity.

     @return CoM velocity (null-angular followed by the linear velocity).
     */
    Eigen::Twistd getVelocity() const;
      
    /**
     Retrieves the acceleration of the CoM.
     
     @return CoM acceleration.
     @warning Currently null.
     */
    Eigen::Twistd getAcceleration() const;
      
    /**
     Retrieves the wrench at the CoM.

     @return Wrench at the CoM
     @warning Currently null.
     */
    Eigen::Wrenchd getWrench() const;
      
      
    /**
     Retrieves the 6xDoF jacobian of the CoM, where the angular part is actually
     null. The linear part is retrieved from the model referenced by this class.

     @return 6-dimensional CoM Jacobian.
     */
    Eigen::Matrix<double,6,Eigen::Dynamic> getJacobian() const;
      
      
    /**
     Indicates whether the frame depends on the mannequein's configuration.

     @return True when it does, false otherwise.
     */
    bool dependsOnModelConfiguration() const;
      
    /**
     Retrieves reference to the model object stored by this class.

     @return Reference to model used to build this object.
     */
    const Model& getModel() const;

      // KDL migration
      /**
       Returns the position of the CoM as stored in the model object passed to this class.
       The orientation is set to the identity.
       
       @return Position of the CoM as a Displacementd object.
       */
      KDL::Frame getPositionKDL() const;
      
      /**
       Retrieves the velocity of the CoM, composed by the linear velocity as stored in the model
       object passed to this class and null angular velocity.
       
       @return CoM velocity (null-angular followed by the linear velocity).
       */
      KDL::Twist getVelocityKDL() const;
      
      /**
       Retrieves the acceleration of the CoM.
       
       @return CoM acceleration.
       @warning Currently null.
       */
      KDL::Twist getAccelerationKDL() const;
      

      /**
       Retrieves the wrench at the CoM.

       @return Wrench at CoM.
       @note Currently nulreturning null wrench.
       */
      KDL::Wrench getWrenchKDL() const;

  private:
    struct Pimpl;
    boost::shared_ptr<Pimpl> pimpl;
  };
}

#endif

// cmake:sourcegroup=Api
