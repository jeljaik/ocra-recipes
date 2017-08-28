#include "ocra/control/Model.h"


namespace ocra
{

  Model::Model(const std::string& name, int ndofs, bool freeRoot, const std::string& jointTorqueVariableName,
               const std::string& forceVariableName, const std::string& configurationVariableName,
               const std::string& internalDofsSuffix, const std::string& externalDofsSuffix)
    : NamedInstance(name)
    , _dofs(ndofs)
    , _fixedRoot(!freeRoot)
  {
    _modelContacts = new ModelContacts(*this);
    if (freeRoot)
    {
      int n_int = _dofs-FREE_DOFS;
      _tau = new BaseVariable(jointTorqueVariableName, n_int);
      BaseVariable* q_root = new BaseVariable(configurationVariableName+externalDofsSuffix, FREE_DOFS);
      BaseVariable* q_int = new BaseVariable(configurationVariableName+internalDofsSuffix, n_int);
      _q = new CompositeVariable(configurationVariableName, *q_root, *q_int);
    }
    else
    {
      _tau = new BaseVariable(jointTorqueVariableName, _dofs);
      _q = new BaseVariable(configurationVariableName, _dofs);
    }
    _q_dot = &_q->getTimeDerivative();
    _q_ddot = &_q_dot->getTimeDerivative();
    _jointDamping = VectorXd::Zero(_tau->getSize());


    _q->connect<EVT_CHANGE_VALUE>(*this, &Model::invalidate);
    _q_dot->connect<EVT_CHANGE_VALUE>(*this, &Model::invalidate);

  }

  Model::~Model()
  {
    _q->disconnect<EVT_CHANGE_VALUE>(*this, &Model::invalidate);
    _q_dot->disconnect<EVT_CHANGE_VALUE>(*this, &Model::invalidate);

    if (!_fixedRoot)
    {
      delete &(*static_cast<CompositeVariable*>(_q))(1);
      delete &(*static_cast<CompositeVariable*>(_q))(0);
    }

    delete _q;
    delete _tau;
    delete _modelContacts;
  }


  int Model::nbDofs() const
  {
    return _dofs;
  }

  int Model::nbInternalDofs() const
  {
    if (!_fixedRoot)
      return _dofs-FREE_DOFS;
    else
      return _dofs;
  }

  bool Model::hasFixedRoot() const
  {
    return _fixedRoot;
  }

  void Model::setJointPositions(const VectorXd& q)
  {
    this->modelMutex.wait();
    getInternalConfigurationVariable().setValue(q);
    doSetJointPositions(q);
    this->modelMutex.post();
  }

  void Model::setJointVelocities(const VectorXd& q_dot)
  {
    this->modelMutex.wait();
    getInternalVelocityVariable().setValue(q_dot);
    doSetJointVelocities(q_dot);
    this->modelMutex.post();
  }

  void Model::setFreeFlyerPosition(const Eigen::Displacementd& H_root)
  {
    assert(!_fixedRoot);
    doSetFreeFlyerPosition(H_root);
    //we do not set H_root in q_root on purpose: H_root is in SE3 while q_root is in R^6
    propagate<EVT_CHANGE_VALUE>();
  }

  void Model::setFreeFlyerVelocity(const Eigen::Twistd& T_root)
  {
    assert(!_fixedRoot);
    doSetFreeFlyerVelocity(T_root);
    getRootVelocityVariable().setValue(T_root.get());
    propagate<EVT_CHANGE_VALUE>();
  }

  void Model::setFreeFlyerPositionKDL(const KDL::Frame& H_root)
  {
      assert(!_fixedRoot);
      doSetFreeFlyerPositionKDL(H_root);
      //we do not set H_root in q_root on purpose: H_root is in SE3 while q_root is in R^6
      propagate<EVT_CHANGE_VALUE>();
  }

  void Model::setFreeFlyerVelocityKDL(const KDL::Twist& T_root)
  {
      assert(!_fixedRoot);
      doSetFreeFlyerVelocityKDL(T_root);
      getRootVelocityVariable().setValue(ocra::util::KDLTwistToEigenVectorXd(T_root));
      propagate<EVT_CHANGE_VALUE>();
  }

  void Model::setState(const VectorXd& q, const VectorXd& q_dot)
  {
    setJointPositions(q);
    setJointVelocities(q_dot);
    doSetState(q, q_dot);
  }

  void Model::setState(const Eigen::Displacementd& H_root, const VectorXd& q, const Eigen::Twistd& T_root, const VectorXd& q_dot)
  {
    setJointPositions(q);
    setJointVelocities(q_dot);
    if(!_fixedRoot)
    {
      setFreeFlyerPosition(H_root);
      setFreeFlyerVelocity(T_root);
    }
    doSetState(H_root, q, T_root, q_dot);
  }

  void Model::setStateKDL(const KDL::Frame& H_root, const VectorXd& q, const KDL::Twist& T_root, const VectorXd& q_dot)
  {
      setJointPositions(q);
      setJointVelocities(q_dot);
      if(!_fixedRoot) // floating base
      {
          setFreeFlyerPositionKDL(H_root);
          setFreeFlyerVelocityKDL(T_root);
      }
      doSetStateKDL(H_root, q, T_root, q_dot);
  }

  Variable& Model::getConfigurationVariable() const
  {
    return *_q;
  }

  Variable& Model::getVelocityVariable() const
  {
    return *_q_dot;
  }

  Variable& Model::getAccelerationVariable() const
  {
    return *_q_ddot;
  }

  Variable& Model::getJointTorqueVariable() const
  {
    return *_tau;
  }

  Variable& Model::getRootConfigurationVariable() const
  {
    assert(!_fixedRoot);
    return (*static_cast<CompositeVariable*>(_q))(0);
  }


  /**
   Returns the internal joint configuration (without floating-base).
   */
  Variable& Model::getInternalConfigurationVariable() const
  {
    if (_fixedRoot)
      return *_q;
    else
      return (*static_cast<CompositeVariable*>(_q))(1);
  }


  /**
   Returns the root velocity variable, which is extracted from the composite
   generalized velocities variable */
  Variable& Model::getRootVelocityVariable() const
  {
    assert(!_fixedRoot);
    return (*static_cast<CompositeVariable*>(_q_dot))(0);
  }


  /**
   Returns the internal joint velocities (without floating-base).
   */
  Variable& Model::getInternalVelocityVariable() const
  {
    if (_fixedRoot)
      return *_q_dot;
    else
      return (*static_cast<CompositeVariable*>(_q_dot))(1);
  }


  /**
   Returns the root acceleration variable, which is extracted from the composite
   generalized accelerations variable.
   */
  Variable& Model::getRootAccelerationVariable() const
  {
    assert(!_fixedRoot);
    return (*static_cast<CompositeVariable*>(_q_ddot))(0);
  }


  /**
   Returns the internal acceleration variable, which is extracted from the composite
   generalized acceleration variable.
   */
  Variable& Model::getInternalAccelerationVariable() const
  {
    if (_fixedRoot)
      return *_q_ddot;
    else
      return (*static_cast<CompositeVariable*>(_q_ddot))(1);
  }


  ModelContacts& Model::getModelContacts() const
  {
    return* _modelContacts;
  }


  /**
   Gets the index corresponding to the name of a given segment.

   @param name Name of the segment.
   @return Integer index of the corresponding segment.
   */
  int Model::getSegmentIndex(const std::string& name) const
  {
    return doGetSegmentIndex(name);
  }


  /**
   Given the integer index of a segment, retrieves the corresponding name.

   @param index Segment index.
   @return Name of the segment.
   */
  const std::string& Model::getSegmentName(int index) const
  {
    assert(index>=0 && index <nbSegments());
    return doGetSegmentName(index);
  }

  void Model::setJointDamping(const VectorXd& damping)
  {
    _jointDamping = damping;
  }

  const VectorXd& Model::getJointDamping() const
  {
    return _jointDamping;
  }

  int Model::doGetSegmentIndex(const std::string& name) const
  {
    throw std::runtime_error("[Model::doGetSegmentIndex] This function was not overriden for a specific model");
  }

  const std::string& Model::doGetSegmentName(int index) const
  {
    throw std::runtime_error("[Model::doGetSegmentName] This function was not overriden for a specific model");
  }

  void Model::invalidate(int timestamp)
  {
    doInvalidate();
    doSetJointPositions(getInternalConfigurationVariable().getValue());
    doSetJointVelocities(getInternalVelocityVariable().getValue());
    if (!_fixedRoot)
      doSetFreeFlyerVelocity(Twistd(getRootVelocityVariable().getValue()));
  }

  void Model::invalidateKDL(int timestamp)
  {
    doInvalidate();
    doSetJointPositions(getInternalConfigurationVariable().getValue());
    doSetJointVelocities(getInternalVelocityVariable().getValue());
    
    //TODO: Create a KDL::Twist from an Eigen::VectorXd (tmpRootVel)
    Eigen::VectorXd tmpRootVel(getRootVelocityVariable().getValue());
    if (!_fixedRoot)
        doSetFreeFlyerVelocityKDL(ocra::util::EigenVectorToKDLTwist(tmpRootVel));
  }

  //TODO: Clean this shit up...
  //-----------------------------BEGINING OF SHIT-------------------------------//
  int Model::getDofIndex(const std::string& name) const
  {
    return doGetDofIndex(name);
  }

  const std::string& Model::getDofName(int index) const
  {
    assert(index>=0 && index < nbDofs());
    return doGetDofName(index);
  }

  const std::string Model::SegmentName(const std::string& name) const
  {
      return doSegmentName(name);
  }

  const std::string Model::DofName(const std::string& name) const
  {
      return doDofName(name);
  }

  int Model::doGetDofIndex(const std::string& name) const
  {
    throw std::runtime_error("[Model::doGetDofIndex] This function was not overriden for a specific model");
  }

  const std::string& Model::doGetDofName(int index) const
  {
    throw std::runtime_error("[Model::doGetDofName] This function was not overriden for a specific model");
  }

  // Translates segment name to the format of the model (for example orcXdeModel should be robot.name)
  const std::string Model::doSegmentName(const std::string& name) const
  {
    throw std::runtime_error("[Model::doSegmentName] This function was not overriden for a specific model");
  }

  // Translates DOF name to the format of the model (for example orcXdeModel should be robot.name)
  const std::string Model::doDofName(const std::string& name) const
  {
    throw std::runtime_error("[Model::doDofName] This function was not overriden for a specific model");
  }
  //-----------------------------END OF SHIT------------------------------------//


}

// cmake:sourcegroup=Api
