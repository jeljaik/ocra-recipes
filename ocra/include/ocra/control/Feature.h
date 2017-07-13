/*!
 \file Feature.h
 \brief A class hierarchy to compute task errors based on control frames.

 Copyright (C) 2010 CEA/DRT/LIST/DTSI/SRCI

 \author Escande Adrien
 \author Evrard Paul
 \date 2010/10/03

 File history:
 Modified in 2017 by Ryan Lober and Jorhabib Eljaik to include KDL and get rid of Eigen LGSM.
 */

#ifndef _OCRA_FEATURE_H_
#define _OCRA_FEATURE_H_

#include "ocra/control/ControlEnum.h"
#include "ocra/optim/NamedInstance.h"

#include "ocra/control/ControlFrame.h"
#include "ocra/control/Model.h"
#include "ocra/control/FullState.h"
#include "ocra/control/PartialState.h"
#include "ocra/control/TaskState.h"
#include <ocra/util/Macros.h>
#include "ocra/util/KDLUtilities.h"

#include <Eigen/Lgsm>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

namespace ocra
{
    class ControlFrame;
    class FullState;
    class PartialState;
}

namespace ocra
{
    // --- ABSTRACT -----------------------------------------------

    /**
     * @class Feature
     * @brief Interface used by tasks to compute errors and jacobians.
     * @details A feature is associated with a control frame and computes errors and jacobians
     *          related to the frame, or to the 'difference' between a frame and another
     *          associated to a different feature, denoted 'desired feature'. A task can be
     *          defined as the error between a feature and a desired feature. Therefore, they
     *          can be seen as key objects to build instances of the Task class.
     * @note See the children class to know the different tasks that can be built.
     * @warning Concurrent calls of computeXXX() and computeXXX(FeatureDes) are not thread-safe
     *          by default (most implementeations use the same cache).
     */
    class Feature : public NamedInstance, boost::noncopyable
    {
        DEFINE_CLASS_POINTER_TYPEDEFS(Feature)
    protected:
        Feature(const std::string& name);

    public:
        virtual ~Feature() = 0;

        virtual int getDimension() const = 0;

        virtual const Eigen::MatrixXd& getSpaceTransform() const = 0;

        virtual const Eigen::VectorXd& computeEffort(const Feature& featureDes) const = 0;
        virtual const Eigen::VectorXd& computeAcceleration(const Feature& featureDes) const = 0;
        virtual const Eigen::VectorXd& computeError(const Feature& featureDes) const = 0;
        virtual const Eigen::VectorXd& computeErrorDot(const Feature& featureDes) const = 0;
        virtual const Eigen::MatrixXd& computeJacobian(const Feature& featureDes) const = 0;
        virtual const Eigen::MatrixXd& computeProjectedMass(const Feature& featureDes) const = 0;
        virtual const Eigen::MatrixXd& computeProjectedMassInverse(const Feature& featureDes) const = 0;

        virtual const Eigen::VectorXd& computeEffort() const = 0;
        virtual const Eigen::VectorXd& computeAcceleration() const = 0;
        virtual const Eigen::VectorXd& computeError() const = 0;
        virtual const Eigen::VectorXd& computeErrorDot() const = 0;
        virtual const Eigen::MatrixXd& computeJacobian() const = 0;
        virtual const Eigen::MatrixXd& computeProjectedMass() const = 0;
        virtual const Eigen::MatrixXd& computeProjectedMassInverse() const = 0;

        virtual TaskState getState() const = 0;
        virtual void setState(const TaskState& newState) = 0;
    };


    // --- POSITION -----------------------------------------------

    /**
     * @class PositionFeature
     * @brief Used to build position tasks.
     * @details The error between two PositionFeature is equal to the difference between
     * the origins of their associated frames (the orientation error is neglected).
     * You can specify axes at construction to consider only some axes: position along
     * x, y, or z; position in the xy, xz, yz planes, or in the cartesian space (x, y, and z).
     * The output space of the jacobian given by this feature correspond to the world frame.
     * Use getSpaceTransform() to get the transformation from the control frame space to
     * the output space.
     */
    class PositionFeature : public Feature
    {
        DEFINE_CLASS_POINTER_TYPEDEFS(PositionFeature)

    public:
        PositionFeature(const std::string& name, ControlFrame::Ptr frame, ECartesianDof axes);

        int getDimension() const;


        /**
         Retrieves the transformation from the control frame space to the output space.

         @return Transformation matrix.
         */
        const Eigen::MatrixXd& getSpaceTransform() const;

        /**
         Difference between the actual force acting on this frame and the desired.

         @param featureDes Desired state of the feature.
         @return 3-dimensional effort.
         */
        const Eigen::VectorXd& computeEffort(const Feature& featureDes) const;
        
        /**
         Given a desired position feature, this method computes the difference between
         the actual frame's linear acceleration and the desired frame's linear acceleration,
         projected on the controlled axes.

         @param featureDes Desired position feature.
         @return Acceleration error in the desired axes.
         */
        const Eigen::VectorXd& computeAcceleration(const Feature& featureDes) const;
        
        /**
         The error of a position feature, given a desired position, is simply the difference between
         the linear positions of the current and desired frames.

         @param featureDes Desired position feature.
         @return Position feature error.
         */
        const Eigen::VectorXd& computeError(const Feature& featureDes) const;
        
        /**
         Computes the difference between the current frame's linear velocity and the desired linear
         velocity of the frame, projected on the controlled axes.

         @param featureDes Desired position feature.
         @return Velocity error of the frame.
         */
        const Eigen::VectorXd& computeErrorDot(const Feature& featureDes) const;
        
        const Eigen::MatrixXd& computeJacobian(const Feature& featureDes) const;
        
        const Eigen::MatrixXd& computeProjectedMass(const Feature& featureDes) const;
        
        const Eigen::MatrixXd& computeProjectedMassInverse(const Feature& featureDes) const;
        
        /**
         Pretty much retrieves the force at the origin of the control frame projected
         on the controlled axes.

         @return 3-dimensional effort.
         */
        const Eigen::VectorXd& computeEffort() const;
        
        /**
         Retrieves the linear acceleration of the current frame projected on the controlled axes.

         @return Frame's linear acceleration.
         */
        const Eigen::VectorXd& computeAcceleration() const;
        
        /**
         Linear position of the current frame projected on the controlled axes.
         @return Frame's linear position.
         */
        const Eigen::VectorXd& computeError() const;
        
        /**
         Linear velocity of the current frame projected on the controlled axes.

         @return Linear velocity error.
         */
        const Eigen::VectorXd& computeErrorDot() const;
        
        const Eigen::MatrixXd& computeJacobian() const;
        
        const Eigen::MatrixXd& computeProjectedMass() const;
        
        const Eigen::MatrixXd& computeProjectedMassInverse() const;
        
        /**
         Returns current state of the frame, composed by its position and velocity.

         @return Current frame's state.
         */
        TaskState getState() const;
        
        /**
         Sets a new state to this position feature.

         @param newState Desired new state.
         */
        void setState(const TaskState& newState);

    private:
        struct Pimpl;
        boost::shared_ptr<Pimpl> pimpl;
    };


    // --- POINT CONTACT ------------------------------------------

    //! Used to build point contact tasks
    /*!
     Desired features have no meaning here, so the versions of the methods
     with desired features passed as arguments will throw. These features
     are used to build point contact tasks, and will typically be associated
     with friction cones. The space transform is identity: the error is
     computed in the controlled frame.
     */
    class PointContactFeature : public Feature
    {
        DEFINE_CLASS_POINTER_TYPEDEFS(PointContactFeature)

    public:
        PointContactFeature(const std::string& name, ControlFrame::Ptr frame);

        int getDimension() const;

        const Eigen::MatrixXd& getSpaceTransform() const;

        const Eigen::VectorXd& computeEffort(const Feature& featureDes) const;
        const Eigen::VectorXd& computeAcceleration(const Feature& featureDes) const;
        const Eigen::VectorXd& computeError(const Feature& featureDes) const;
        const Eigen::VectorXd& computeErrorDot(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeJacobian(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeProjectedMass(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeProjectedMassInverse(const Feature& featureDes) const;

        const Eigen::VectorXd& computeEffort() const;
        const Eigen::VectorXd& computeAcceleration() const;
        const Eigen::VectorXd& computeError() const;
        const Eigen::VectorXd& computeErrorDot() const;
        const Eigen::MatrixXd& computeJacobian() const;
        const Eigen::MatrixXd& computeProjectedMass() const;
        const Eigen::MatrixXd& computeProjectedMassInverse() const;

        TaskState getState() const;
        void setState(const TaskState& newState);

    private:
        struct Pimpl;
        boost::shared_ptr<Pimpl> pimpl;
    };


    // --- ORIENTATION --------------------------------------------

    /**
     * @class OrientationFeature
     * @brief Used to build orientation tasks.
     * @details OrientationFeature instances consider only the orientation of the associated control frame.
     * The error between two orientation features is the log of Rdes.inverse() * R, where R is
     * the orientation of the orientation frame and Rdes is the desired orientation.
     */
    class OrientationFeature : public Feature
    {
        DEFINE_CLASS_POINTER_TYPEDEFS(OrientationFeature)

    public:
        OrientationFeature(const std::string& name, ControlFrame::Ptr frame);

        int getDimension() const;

        const Eigen::MatrixXd& getSpaceTransform() const;
        /**
         Given a desired orientation feature, it computes the effort as the difference between the torque at the frame for the current orientation and the one needed for the desired orientation.

         @param featureDes Desired orientation feature.
         @return 3-dimensional effort on each axis of rotation.
         */
        const Eigen::VectorXd& computeEffort(const Feature& featureDes) const;

        /**
         Given a desired orientation feature, it computes the angular acceleration that the associated control frame need to go from the current orientation to the desired one.

         @param featureDes Desired orientation feature.
         @return Angular acceleration needed.
         */
        const Eigen::VectorXd& computeAcceleration(const Feature& featureDes) const;

        /**
         Given a desired orientation feature, returns the logarithmic error in the orientation.

         @param featureDes Desired orientation feature.
         @return Orientation error (logarithm)
         */
        const Eigen::VectorXd& computeError(const Feature& featureDes) const;

        /**
         Computes the angular velocity needed to go from the current configuration to the desired one.

         @param featureDes Desired orientation feature.
         @return Angular velocity error.
         */
        const Eigen::VectorXd& computeErrorDot(const Feature& featureDes) const;

        const Eigen::MatrixXd& computeJacobian(const Feature& featureDes) const;

        const Eigen::MatrixXd& computeProjectedMass(const Feature& featureDes) const;

        const Eigen::MatrixXd& computeProjectedMassInverse(const Feature& featureDes) const;

        /**
         Returns the torque at the current control frame.

         @return Current effort.
         */
        const Eigen::VectorXd& computeEffort() const;

        /**
         Pretty much the current angular acceleration associated to this control frame.

         @return Angular acceleration of the control frame
         */
        const Eigen::VectorXd& computeAcceleration() const;

        /**
         Returns the logarithmic error in the orientation

         @return Orientation error.
         */
        const Eigen::VectorXd& computeError() const;

        /**
         Pretty much returns the angular velocity of the current control frame.

         @return Angular velocity.
         */
        const Eigen::VectorXd& computeErrorDot() const;

        const Eigen::MatrixXd& computeJacobian() const;

        const Eigen::MatrixXd& computeProjectedMass() const;

        const Eigen::MatrixXd& computeProjectedMassInverse() const;
        /**
         Retrieves the current state of the associated control frame. The state is composed by its position, velocity, acceleration and wrench.

         @return State of the feature.
         */
        TaskState getState() const;

        /**
         Sets the state of the feature, i.e. position, velocity, acceleration and wrench.

         @param newState New state to be set.
         */
        void setState(const TaskState& newState);

    private:
        struct Pimpl;
        boost::shared_ptr<Pimpl> pimpl;
    };


    // --- DISPLACEMENT -------------------------------------------

    //! Used to build position/orientation tasks.
    /*!
     These features combine position and orientation features, to control both the position
     and orientation of a control frame. See the documentation of PositionFeature and
     OrientationFeature.
     */
    class DisplacementFeature : public Feature
    {
        DEFINE_CLASS_POINTER_TYPEDEFS(DisplacementFeature)

    public:
        DisplacementFeature(const std::string& name, ControlFrame::Ptr frame, ECartesianDof axes);

        int getDimension() const;

        const Eigen::MatrixXd& getSpaceTransform() const;

        const Eigen::VectorXd& computeEffort(const Feature& featureDes) const;
        const Eigen::VectorXd& computeAcceleration(const Feature& featureDes) const;
        const Eigen::VectorXd& computeError(const Feature& featureDes) const;
        const Eigen::VectorXd& computeErrorDot(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeJacobian(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeProjectedMass(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeProjectedMassInverse(const Feature& featureDes) const;

        const Eigen::VectorXd& computeEffort() const;
        const Eigen::VectorXd& computeAcceleration() const;
        const Eigen::VectorXd& computeError() const;
        const Eigen::VectorXd& computeErrorDot() const;
        const Eigen::MatrixXd& computeJacobian() const;
        const Eigen::MatrixXd& computeProjectedMass() const;
        const Eigen::MatrixXd& computeProjectedMassInverse() const;

        TaskState getState() const;
        void setState(const TaskState& newState);

    private:
        struct Pimpl;
        boost::shared_ptr<Pimpl> pimpl;
    };


    // --- CONTACT CONSTRAINT -------------------------------------

    //! Used to build contact constraint tasks (null acceleration).
    /*!
     These features is similar to a displacement feature, except all quantities are
     expressed in the controlled frame (space transform is identity). Desired
     features are irrelevant and computeXXX(const Feature&) will throw.
     */
    class ContactConstraintFeature : public Feature
    {
        DEFINE_CLASS_POINTER_TYPEDEFS(ContactConstraintFeature)

    public:
        ContactConstraintFeature(const std::string& name, ControlFrame::Ptr frame);

        int getDimension() const;

        const Eigen::MatrixXd& getSpaceTransform() const;

        const Eigen::VectorXd& computeEffort(const Feature& featureDes) const;
        const Eigen::VectorXd& computeAcceleration(const Feature& featureDes) const;
        const Eigen::VectorXd& computeError(const Feature& featureDes) const;
        const Eigen::VectorXd& computeErrorDot(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeJacobian(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeProjectedMass(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeProjectedMassInverse(const Feature& featureDes) const;

        const Eigen::VectorXd& computeEffort() const;
        const Eigen::VectorXd& computeAcceleration() const;
        const Eigen::VectorXd& computeError() const;
        const Eigen::VectorXd& computeErrorDot() const;
        const Eigen::MatrixXd& computeJacobian() const;
        const Eigen::MatrixXd& computeProjectedMass() const;
        const Eigen::MatrixXd& computeProjectedMassInverse() const;


        TaskState getState() const;
        void setState(const TaskState& newState);


    private:
        struct Pimpl;
        boost::shared_ptr<Pimpl> pimpl;
    };





    // --- ARTICULAR ----------------------------------------------

    //! Used to build tasks in the manikin configuration space.
    /*!
     Can be used e.g. for Joint PD control.
     */
    class FullStateFeature : public Feature
    {
        DEFINE_CLASS_POINTER_TYPEDEFS(FullStateFeature)

    public:
        FullStateFeature(const std::string& name, FullState::Ptr state);

        int getDimension() const;

        const Eigen::MatrixXd& getSpaceTransform() const;

        const Eigen::VectorXd& computeEffort(const Feature& featureDes) const;
        const Eigen::VectorXd& computeAcceleration(const Feature& featureDes) const;
        const Eigen::VectorXd& computeError(const Feature& featureDes) const;
        const Eigen::VectorXd& computeErrorDot(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeJacobian(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeProjectedMass(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeProjectedMassInverse(const Feature& featureDes) const;

        const Eigen::VectorXd& computeEffort() const;
        const Eigen::VectorXd& computeAcceleration() const;
        const Eigen::VectorXd& computeError() const;
        const Eigen::VectorXd& computeErrorDot() const;
        const Eigen::MatrixXd& computeJacobian() const;
        const Eigen::MatrixXd& computeProjectedMass() const;
        const Eigen::MatrixXd& computeProjectedMassInverse() const;


        TaskState getState() const;
        void setState(const TaskState& newState);


    private:
        struct Pimpl;
        boost::shared_ptr<Pimpl> pimpl;
    };

    // --- PARTIAL - ARTICULAR ----------------------------------------

    //! Used to build tasks in a partial configuration space.
    /*!
     Can be used e.g. for Joint PD control.
     */
    class PartialStateFeature : public Feature
    {
        DEFINE_CLASS_POINTER_TYPEDEFS(PartialStateFeature)

    public:
        PartialStateFeature(const std::string& name, PartialState::Ptr state);

        int getDimension() const;

        const Eigen::MatrixXd& getSpaceTransform() const;

        const Eigen::VectorXd& computeEffort(const Feature& featureDes) const;
        const Eigen::VectorXd& computeAcceleration(const Feature& featureDes) const;
        const Eigen::VectorXd& computeError(const Feature& featureDes) const;
        const Eigen::VectorXd& computeErrorDot(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeJacobian(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeProjectedMass(const Feature& featureDes) const;
        const Eigen::MatrixXd& computeProjectedMassInverse(const Feature& featureDes) const;

        const Eigen::VectorXd& computeEffort() const;
        const Eigen::VectorXd& computeAcceleration() const;
        const Eigen::VectorXd& computeError() const;
        const Eigen::VectorXd& computeErrorDot() const;
        const Eigen::MatrixXd& computeJacobian() const;
        const Eigen::MatrixXd& computeProjectedMass() const;
        const Eigen::MatrixXd& computeProjectedMassInverse() const;


        TaskState getState() const;
        void setState(const TaskState& newState);


    private:
        struct Pimpl;
        boost::shared_ptr<Pimpl> pimpl;
    };
}

#endif

// cmake:sourcegroup=Api
