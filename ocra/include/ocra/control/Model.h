/** @file Model.h
  * @brief Declaration file of the Model class.
  *
  *   Copyright (C) 2010 CEA/DRT/LIST/DIASI/LSI
  *
  * @author Escande Adrien
  *	@date 10/08/02
  */

#ifndef _OCRACONTROL_MODEL_H_
#define _OCRACONTROL_MODEL_H_

// includes
#include <ocra/util/MathTypes.h>
#include "ocra/optim/ObserverSubject.h"
#include "ocra/optim/Variable.h"
#include "ocra/optim/NamedInstance.h"
#include "ocra/control/ModelContacts.h"
#include <ocra/util/Macros.h>
#include <ocra/util/KDLUtilities.h>
#include <Eigen/Lgsm>
#include "kdl/frames_io.hpp"
#include "kdl/frames.hpp"


#include <string>
#include <iostream>

#include <yarp/os/Semaphore.h>

namespace ocra
{
  /** @class Model
    *	@brief %Model class.
    *	@warning None
    *
    *   @param name Name of the model.
    *   @param ndofs Total number of degrees of freedom of the model. Includes floating-base.
    *   @param freeRoot FALSE for floating-base robots. True otherwise.
    *   @param jointTorqueVariableName Name given to the joint torques variable.
    *   @param forceVariableName Name given to the force variables.
    *   @param configurationVariableName Name given to the joint configuration (angles) variable.
    *   @param internalDofsSuffix Name Suffix appended to the name of the internal degrees of freedom variable (no floating-base).
    *   @param externalDofsSuffix Name Suffix appended to the name of the floating-base DOF variable.
    *
    *
    * TODO: complete description
    *
    * terms of the dynamic equation are given so that the equation writes this way :
    * \f$ M\dot{T} + N + G = L tau - J_c^T f \f$
    */
  class Model : public ObserverSubject, public NamedInstance
  {
      DEFINE_CLASS_POINTER_TYPEDEFS(Model)
  public:
    Model(const std::string& name, int ndofs, bool freeRoot,
          const std::string& jointTorqueVariableName = "tau",
          const std::string& forceVariableName = "f",
          const std::string& configurationVariableName = "q",
          const std::string& internalDofsSuffix = "_int",
          const std::string& externalDofsSuffix = "_root");
  public:
    virtual ~Model();

    // ------------------------ public interface --------------------------------
  public:
    //set/get configuration
    void  setJointPositions(const Eigen::VectorXd& q);
    void  setJointVelocities(const Eigen::VectorXd& q_dot);
    void  setFreeFlyerPosition(const Eigen::Displacementd& H_root);
    void  setFreeFlyerVelocity(const Eigen::Twistd& T_root);
    void  setFreeFlyerPositionKDL(const KDL::Frame& H_root);
    void  setFreeFlyerVelocityKDL(const KDL::Twist& T_root);
    void  setState(const Eigen::VectorXd& q, const Eigen::VectorXd& q_dot);
    void  setState(const Eigen::Displacementd& H_root, const Eigen::VectorXd& q, const Eigen::Twistd& T_root, const Eigen::VectorXd& q_dot);
    void  setStateKDL(const KDL::Frame& H_root, const Eigen::VectorXd& q, const KDL::Twist& T_root, const Eigen::VectorXd& q_dot);
    virtual const Eigen::VectorXd& getJointPositions()         const = 0;
    virtual const Eigen::VectorXd& getJointVelocities()        const = 0;
    virtual const Eigen::VectorXd& getJointAccelerations()     const { std::cout << "getJointAccelerations() has not been implemented" << std::endl; };
    virtual const Eigen::VectorXd& getJointTorques()           const = 0;
    virtual const Eigen::Displacementd& getFreeFlyerPosition() const = 0;
    virtual const Eigen::Twistd& getFreeFlyerVelocity()        const = 0;
    virtual const KDL::Frame& getFreeFlyerPositionKDL() const = 0;
    virtual const KDL::Twist& getFreeFlyerVelocityKDL()        const = 0;

    //get whole body data
      //dofs
    int                       nbDofs()              const;
    /** Number of internal degrees of freedom, i.e. without the 6 DOF of the floating base*/
    int                       nbInternalDofs()      const;
    bool                      hasFixedRoot()        const;
    virtual int               nbSegments()          const = 0;
    virtual const Eigen::VectorXd&   getActuatedDofs()     const = 0;
    virtual const Eigen::VectorXd&   getJointLowerLimits() const = 0;
    virtual const Eigen::VectorXd&   getJointUpperLimits() const = 0;
      //CoM
    virtual double            getMass()             const = 0;
    virtual const Eigen::Vector3d&   getCoMPosition()      const = 0;
    virtual const Eigen::Vector3d&   getCoMVelocity()      const = 0;
    virtual const Eigen::Vector3d&   getCoMAcceleration()  const { std::cout << "getCoMAcceleration() Not implemented" << std::endl; };
    virtual const Eigen::Vector3d&   getCoMAngularVelocity()      const{ std::cout << "getCoMAngularVelocity() Not implemented" << std::endl; }
    virtual const Eigen::Vector3d&   getCoMJdotQdot()      const = 0;
    virtual const Eigen::Matrix<double,3,Eigen::Dynamic>& getCoMJacobian()      const = 0;
    virtual const Eigen::Matrix<double,3,Eigen::Dynamic>& getCoMAngularJacobian() const{ std::cout << "getCoMAngularVelocity() Not implemented" << std::endl; }
    virtual const Eigen::Matrix<double,3,Eigen::Dynamic>& getCoMJacobianDot()   const = 0;
      //dynamic/static equation terms
    virtual const Eigen::MatrixXd&   getInertiaMatrix()            const = 0;
    virtual const Eigen::MatrixXd&   getInertiaMatrixInverse()     const = 0;
    virtual const Eigen::MatrixXd&   getDampingMatrix()            const = 0;
    virtual const Eigen::VectorXd&   getNonLinearTerms()           const = 0;
    virtual const Eigen::VectorXd&   getLinearTerms()              const = 0;
    virtual const Eigen::VectorXd&   getGravityTerms()             const = 0;

    //segment data
    virtual const Eigen::Displacementd&  getSegmentPosition(int index) const = 0;
    virtual const Eigen::Twistd&         getSegmentVelocity(int index) const = 0;
    virtual const KDL::Frame&            getSegmentPositionKDL(int index) const = 0;
    virtual const KDL::Twist&            getSegmentVelocityKDL(int index) const = 0;
    virtual const Eigen::Matrix<double,6,Eigen::Dynamic>&     getSegmentJacobian(int index) const = 0;
    virtual const Eigen::Matrix<double,6,Eigen::Dynamic>&     getSegmentJdot(int index)     const = 0;
    virtual const Eigen::Twistd&         getSegmentJdotQdot(int index) const = 0;
      virtual const KDL::Twist&          getSegmentJdotQdotKDL(int index) const = 0;
    virtual const Eigen::Matrix<double,6,Eigen::Dynamic>& getJointJacobian(int index) const = 0;
    virtual double                       getSegmentMass(int index) const = 0;
    virtual const Eigen::Vector3d&       getSegmentCoM(int index) const = 0;
    virtual const Eigen::Matrix<double, 6, 6>& getSegmentMassMatrix(int index) const = 0;
    virtual const Eigen::Vector3d&       getSegmentMomentsOfInertia(int index) const = 0;
    virtual const Eigen::Rotation3d&     getSegmentInertiaAxes(int index) const = 0;

    //segment data
    const Eigen::Displacementd&  getSegmentPosition(const std::string& segName) const {
        return getSegmentPosition(getSegmentIndex(segName));
    }
    const Eigen::Twistd& getSegmentVelocity(const std::string& segName) const {
        return getSegmentVelocity(getSegmentIndex(segName));
    }
    const KDL::Frame&  getSegmentPositionKDL(const std::string& segName) const {
        return getSegmentPositionKDL(getSegmentIndex(segName));
    }
    const KDL::Twist&  getSegmentVelocityKDL(const std::string& segName) const {
        return getSegmentVelocityKDL(getSegmentIndex(segName));
    }
    const Eigen::Matrix<double,6,Eigen::Dynamic>&     getSegmentJacobian(const std::string& segName) const {
        return getSegmentJacobian(getSegmentIndex(segName));
    }
    const Eigen::Matrix<double,6,Eigen::Dynamic>&     getSegmentJdot(const std::string& segName)     const {
        return getSegmentJdot(getSegmentIndex(segName));
    }
    const Eigen::Twistd&         getSegmentJdotQdot(const std::string& segName) const {
        return getSegmentJdotQdot(getSegmentIndex(segName));
    }
    const KDL::Twist&         getSegmentJdotQdotKDL(const std::string& segName) const {
        return getSegmentJdotQdotKDL(getSegmentIndex(segName));
    }
    const Eigen::Matrix<double,6,Eigen::Dynamic>& getJointJacobian(const std::string& segName) const {
        return getJointJacobian(getSegmentIndex(segName));
    }
    double                       getSegmentMass(const std::string& segName) const {
        return getSegmentMass(getSegmentIndex(segName));
    }
    const Eigen::Vector3d&       getSegmentCoM(const std::string& segName) const {
        return getSegmentCoM(getSegmentIndex(segName));
    }
    const Eigen::Matrix<double, 6, 6>& getSegmentMassMatrix(const std::string& segName) const {
        return getSegmentMassMatrix(getSegmentIndex(segName));
    }
    const Eigen::Vector3d&       getSegmentMomentsOfInertia(const std::string& segName) const {
        return getSegmentMomentsOfInertia(getSegmentIndex(segName));
    }
    const Eigen::Rotation3d&     getSegmentInertiaAxes(const std::string& segName) const {
        return getSegmentInertiaAxes(getSegmentIndex(segName));
    }

    void setJointDamping(const Eigen::VectorXd& damping);
    const Eigen::VectorXd& getJointDamping() const;

    //variables
    Variable& getConfigurationVariable()  const;
    Variable& getVelocityVariable()       const;
    Variable& getAccelerationVariable()   const;
    Variable& getJointTorqueVariable()    const;

    Variable& getRootConfigurationVariable()      const;
    Variable& getInternalConfigurationVariable()  const;
    Variable& getRootVelocityVariable()           const;
    Variable& getInternalVelocityVariable()       const;
    Variable& getRootAccelerationVariable()       const;
    Variable& getInternalAccelerationVariable()   const;
    yarp::os::Semaphore modelMutex;

    //subModels
    ModelContacts&       getModelContacts() const;

    //utils
    int getSegmentIndex(const std::string& name) const;
    const std::string& getSegmentName(int index) const;

    // ------------------------ public methods ----------------------------------
  public:
      //TODO: Clean this shit up...
      //-----------------------------BEGINING OF SHIT-------------------------------//
      int getDofIndex(const std::string& name) const;
      const std::string& getDofName(int index) const;
      const std::string DofName(const std::string& name) const;
      const std::string SegmentName(const std::string& name) const;
      virtual const std::string&           getJointName             (int index) const = 0;

      //-----------------------------END OF SHIT------------------------------------//


    // ------------------------ public static methods ---------------------------
  public:

    // ------------------------ protected methods -------------------------------
  protected:
    virtual void  doSetState(const Eigen::VectorXd& q, const Eigen::VectorXd& q_dot){};
    virtual void  doSetState(const Eigen::Displacementd& H_root, const Eigen::VectorXd& q, const Eigen::Twistd& T_root, const Eigen::VectorXd& q_dot){};
    virtual void  doSetStateKDL(const KDL::Frame& H_root, const Eigen::VectorXd& q, const KDL::Twist& T_root, const Eigen::VectorXd& q_dot){};

    virtual void  doSetJointPositions(const Eigen::VectorXd& q) = 0;
    virtual void  doSetJointVelocities(const Eigen::VectorXd& q_dot) = 0;
    virtual void  doSetFreeFlyerPosition(const Eigen::Displacementd& H_root) = 0;
    virtual void  doSetFreeFlyerVelocity(const Eigen::Twistd& T_root) = 0;
    virtual void  doSetFreeFlyerPositionKDL(const KDL::Frame& H_root) = 0;
    virtual void  doSetFreeFlyerVelocityKDL(const KDL::Twist& T_root) = 0;

    virtual int                 doGetSegmentIndex(const std::string& name)  const = 0;
    virtual const std::string&  doGetSegmentName(int index)                 const = 0;

    virtual void doInvalidate() {}

    //TODO: Clean this shit up...
    //-----------------------------BEGINING OF SHIT-------------------------------//
    virtual int                 doGetDofIndex       (const std::string& name) const;
    virtual const std::string&  doGetDofName        (int index) const;
    virtual const std::string   doSegmentName       (const std::string& name) const;
    virtual const std::string   doDofName           (const std::string& name) const;
    //-----------------------------END OF SHIT------------------------------------//


    // ------------------------ protected static methods ------------------------
  protected:

    // ------------------------ private methods ---------------------------------
  private:
    void invalidate(int timestamp);
    void invalidateKDL(int timestamp);
    // ------------------------ private static methods --------------------------
  private:

    // ------------------------ protected members -------------------------------
  protected:

    // ------------------------ protected static members ------------------------
  protected:

    // ------------------------ private members ---------------------------------
  private:

    /** Total number of degrees of freedom, floating-base included */
    int         _dofs;
    /** When FALSE, we're dealing with a floating-base robot such as a biped */
    bool        _fixedRoot;
    VectorXd    _jointDamping;
    /** Generalized joint torques variable (floating-base included) */
    Variable*   _tau;
    /** Generalized joint angles variable (floating-base included) */
    Variable*   _q;
    /** Genralized joint velocities variable (floating-base included) */
    Variable*   _q_dot;
    /** Generalized joint accelerations variable (floating-base included) */
    Variable*   _q_ddot;

    ModelContacts* _modelContacts;
    // ------------------------ private static members --------------------------
  private:
    static const int   FREE_DOFS = 6;

    // ------------------------ friendship declarations -------------------------
  };
}


#endif //_OCRACONTROL_MODEL_H_

// cmake:sourcegroup=Api
