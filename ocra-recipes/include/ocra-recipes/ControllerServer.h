/*! \file       ControllerServer.h
 *  \brief      Module base class for the controller Server.
 *  \class      ControllerServer
 *  \details
 *  \author     [Ryan Lober](http://www.ryanlober.com)
 *  \copyright  GNU General Public License.
 */
/*
 *  This file is part of ocra-recipes.
 *  Copyright (C) 2016 Institut des Systèmes Intelligents et de Robotique (ISIR)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */



#ifndef CONTROLLER_SERVER_H
#define CONTROLLER_SERVER_H

#include <iostream>
#include <memory>

#include <ocra/util/Macros.h>
#include <ocra/control/Controller.h>
#include <ocra/control/Model.h>
#include <ocra/control/TaskBuilders/TaskConstructionManager.h>

#include <ocra/optim/OneLevelSolver.h>
#include <wocra/WocraController.h>
#include <hocra/HocraController.h>

#include <Eigen/Dense>
#include <Eigen/Lgsm>
#include <ocra/util/FileOperations.h>

// TODO: Should put in defines for yarp independent builds
#include <yarp/os/Bottle.h>
#include <yarp/os/Port.h>

#include <ocra-recipes/ServerCommunications.h>
#include <ocra-recipes/RobotState.h>

namespace ocra_recipes
{

enum CONTROLLER_TYPE
{
    WOCRA_CONTROLLER = 1,
    HOCRA_CONTROLLER,
    GOCRA_CONTROLLER
};

enum SOLVER_TYPE
{
    QUADPROG = 1,
    QPOASES
};

class ControllerServer
{
    DEFINE_CLASS_POINTER_TYPEDEFS(ControllerServer)
protected:
    
    /**
     The user must implement this method for the robot-specific server. It provides access to the model information of the robot.

     @return Pointer to the model class.
     */
    virtual std::shared_ptr<Model> loadRobotModel() = 0;
    
    /**
     Retrieves the robot state defined by its joints' configuration, joints' velociites, roto-translation matrix of the root reference frame and the twist of the robot root frame. This is specific to every robot. For iCub for example, the WholeBodyInterface class is used as a yarp interface with the robot.

     @param[out] q Joint configuration.
     @param[out] qd Joint velocities.
     @param[out] H_root Root roto-translation matrix.
     @param[out] T_root Root twist vector.
     */
    virtual void getRobotState(Eigen::VectorXd& q, Eigen::VectorXd& qd, Eigen::Displacementd& H_root, Eigen::Twistd& T_root) = 0;
    
    /**
     Retrieves the robot state defined by its joints' configuration, joints' velociites, roto-translation matrix of the root reference frame and the twist of the robot root frame.

     @param[out] q Joint configuration.
     @param[out] qd Joint velocities.
     @param[out] H_root Root roto-translation matrix.
     @param[out] T_root Root twist vector.
     */
    virtual void getRobotStateKDL(Eigen::VectorXd& q, Eigen::VectorXd& qd, KDL::Frame& H_root, KDL::Twist& T_root) = 0;

public:
    /**
     Constructor. 
     
     @param[in] ctrlType Controller type (e.g. WOCRA_CONTROLLER).
     @param[in] solver Solver type (e.g. QUADPROG).
     @param[in] usingInterprocessCommunication [true] Is interprocess communication used. 
     @param[in] useOdometry [false] True if the robot won't remain still during the tasks performed (in walking for example).
     */
    ControllerServer(CONTROLLER_TYPE ctrlType/*=WOCRA_CONTROLLER*/,
                     SOLVER_TYPE solver/*=QUADPROG*/,
                     bool usingInterprocessCommunication=true,
                     bool useOdometry=false);
    /** 
     Default destructor
     */
    virtual ~ControllerServer();

    /**
     Performs initialization tasks (robot-agnostic) in the following order:
       - Load robot model
       - Create object to store robot state.
       - Set the solver.
       - Construct the desired controller (e.g. WOCRA).
       - Activate inter-process communication if specified. 
       - Update model if odometry is not used.

     @return True when all initialization steps are successful.
     */
    bool initialize();
    
    /**
     Actually compute actuation torques.

     @return Vector of internal torques that satisfy the tasks and constraints of the problem.
     */
    const Eigen::VectorXd& computeTorques();
    
    /**
     Updates robot model and then actually computes joint actuation torques.

     @param[out] torques Reference to variable storing computed internal torques.
     */
    void computeTorques(Eigen::VectorXd& torques);

    /**
     Provide access to internal controller.

     @return Pointer to specified controller.
     */
    const std::shared_ptr<ocra::Controller> getController(){return controller;}
    
    /**
     Provide access to internal model.

     @return Pointer to model.
     */
    const std::shared_ptr<ocra::Model> getRobotModel(){return model;}

    bool addTasksFromXmlFile(const std::string& filePath);
    
    bool addTasks(std::vector<ocra::TaskBuilderOptions>& tmOpts);

    /**
     If the useOdometry flag is passed to the server, odometry is computed.

     @return True if properly initialized, false otherwise.
     */
    bool initializeOdometry();
    
    
    /**
     Performs the following steps:
      - Calls the robot-specific method `getRobotState()`.
      - Sets the robot state in the internal model.
      - Broadcasts the robot state through yarp ports.
     */
    void updateModel();
    
    /**
     Performs the following steps:
     - Calls the robot-specific method `getRobotState()`.
     - Sets the robot state in the internal model.
     - Broadcasts the robot state through yarp ports.
     */
    void updateModelKDL();

    /**
     Sets regularization weights to the internal solver.

     @param wDdq
     @param wTau
     @param wFc 
     @todo Ryan needs to give a better description for the input parameters.
     */
    void setRegularizationTermWeights(double wDdq, double wTau, double wFc);

protected:

    ocra::Model::Ptr                     model;
    ocra::Controller::Ptr           controller;
    ocra::Solver::Ptr           internalSolver;

    bool firstRun;

    ServerCommunications::Ptr       serverComs;

    Eigen::VectorXd            tau;
    RobotState              rState;
    Eigen::VectorXd         qdPrevious;
    Eigen::VectorXd         qPrevious;

    CONTROLLER_TYPE    controllerType;
    SOLVER_TYPE            solverType;
    bool                    usingComs;
    bool                  usingOdometry;

    yarp::os::Bottle statesBottle;
    yarp::os::Port statesPort;

    //Odometry options
    std::string modelFile;
    std::string initialFixedFrame;
};

} // namespace ocra_recipes
#endif // CONTROLLER_SERVER_H
