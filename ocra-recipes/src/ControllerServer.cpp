#include <ocra-recipes/ControllerServer.h>

using namespace ocra_recipes;

ControllerServer::ControllerServer(CONTROLLER_TYPE ctrlType,
                                   SOLVER_TYPE solver,
                                   bool usingInterprocessCommunication,
                                   bool useOdometry)
: controllerType(ctrlType)
, solverType(solver)
, usingComs(usingInterprocessCommunication)
, usingOdometry(useOdometry)
{
}

ControllerServer::~ControllerServer()
{
    // if(taskManagerSet)
    //     taskManagerSet->clearSet();
}

bool ControllerServer::initialize()
{
    bool res = true;
    model = loadRobotModel();
    firstRun = true;

    rState = RobotState(model->nbDofs());
    rState.T_root = Eigen::Twistd::Zero();
    if(model)
    {
        // Set the solver.
        switch (solverType)
        {
            case QUADPROG:
            {
                std::cout << "Using QuadProg++ "<<std::endl;
                internalSolver = std::make_shared<ocra::OneLevelSolverWithQuadProg>();
            }break;

            default:
            {
                std::cout << "Using QuadProg++ "<<std::endl;
                internalSolver = std::make_shared<ocra::OneLevelSolverWithQuadProg>();
            }break;
        }

        #ifdef USING_QPOASES
        if(solverType == QPOASES)
        {
            std::cout << "Using qpOASES"<<std::endl;
            internalSolver = std::make_shared<ocra::OneLevelSolverWithQPOASES>();
        }
        #endif //USING_QPOASES


        // Construct the desired controller.
        switch (controllerType)
        {
            case WOCRA_CONTROLLER:
            {
                std::cout << "Constructing a WOCRA Controller "<<std::endl;
                bool useReducedProblem = false;
                controller = std::make_shared<wocra::WocraController>("WocraController", model, std::static_pointer_cast<ocra::OneLevelSolver>(internalSolver), useReducedProblem);
            }break;

            case HOCRA_CONTROLLER:
            {
                std::cout << "Constructing a HOCRA Controller "<<std::endl;
                bool useReducedProblem = false;
                controller = std::make_shared<hocra::HocraController>("HocraController", model, std::static_pointer_cast<ocra::OneLevelSolver>(internalSolver), useReducedProblem);
            }break;

            default:
            {
                std::cout << "Constructing a WOCRA Controller [Default]"<<std::endl;
                bool useReducedProblem = false;
                controller = std::make_shared<wocra::WocraController>("WocraController", model, std::static_pointer_cast<ocra::OneLevelSolver>(internalSolver), useReducedProblem);
            }break;
        }
    }

    if(usingComs)
    {
        serverComs = std::make_shared<ServerCommunications>(controller, model);
        res &= serverComs->open();
        res &= statesPort.open("/ControllerServer/states:o");
    }

    res &= bool(model);
    res &= bool(controller);

    // WARNING! If useOdometry is true we must call updateModel during initialization explicitly after ControllerServer::initialize()
    if (!this->usingOdometry) {
        // Setting up initial contact state. Both feet on the ground.
        // TODO: The initial contact state should be obtained from the initial configuration of the robot automatically.
        if (!model->hasFixedRoot()) {
            this->controller->setContactState(1,1);
        }
#ifndef OCRA_USES_KDL
        updateModel();
#else
        updateModelKDL();
#endif
    }
    return res;
}

const Eigen::VectorXd& ControllerServer::computeTorques()
{
    computeTorques(tau);
    return tau;
}

void ControllerServer::computeTorques(Eigen::VectorXd& torques)
{
#ifndef OCRA_USES_KDL
    updateModel();
#else
    updateModelKDL();
#endif
    controller->computeOutput(torques);
}

void ControllerServer::updateModel()
{
    getRobotState(rState.q, rState.qd, rState.H_root, rState.T_root);
    if (model->hasFixedRoot()){
        model->setState(rState.q, rState.qd);
    }else{
//         std::string home = "/home/jorhabib/Documents/debugging";
//         std::string twistHome = std::string(home + "/fbTwist");
//         std::string positionHome = std::string(home + "/fbPosition");
        model->setState(rState.H_root, rState.q, rState.T_root, rState.qd);
//         ocra::utils::writeToFile(rState.T_root, twistHome);
//         Eigen::Vector3d transTmp = rState.H_root.getTranslation();
//         ocra::utils::writeToFile(transTmp, positionHome);
//         ocra::utils::writeToFile(rState.T_root, twistHome);
    }
    if (!statesPort.write( rState)) {
        OCRA_ERROR("Couldn't write robot state for client. Not really doing anything about it, except reporting.");
    }
}

void ControllerServer::updateModelKDL()
{
    getRobotStateKDL(rState.q, rState.qd, rState.H_rootKDL, rState.T_rootKDL);
    Eigen::VectorXd troot(6);
    troot = ocra::util::KDLTwistToEigenVectorXd(rState.T_rootKDL);
    //TODO: FOR DEBUGGING PURPOSES. REMOVE THESE LINES IF IT'S WRONG.
    Eigen::VectorXd tmpRoot(6); tmpRoot.setZero();
    rState.T_rootKDL = KDL::Twist::Zero();
    Eigen::VectorXd tmpRootPos(3);
    tmpRootPos = ocra::util::KDLVectorToEigenVector3d(rState.H_rootKDL.p);
    std::string home = "/home/jorhabib/Documents/debugging";
    std::string twistHome = std::string(home + "/fbTwist");
    std::string positionHome = std::string(home + "/fbPosition");
    ocra::utils::writeToFile(troot,twistHome);
    ocra::utils::writeToFile(tmpRootPos,positionHome);
    if (model->hasFixedRoot()){
        model->setState(rState.q, rState.qd);
    } else {
        model->setStateKDL(rState.H_rootKDL, rState.q, rState.T_rootKDL, rState.qd);
    }
    if (!statesPort.write(rState)) {
        OCRA_ERROR("Couldn't write robot state for client. Not really doing anything about it, except reporting.");
    }
}

bool ControllerServer::addTasksFromXmlFile(const std::string& filePath)
{
    ocra::TaskConstructionManager factory(model, controller, filePath);
    return true;
}

bool ControllerServer::addTasks(std::vector<ocra::TaskBuilderOptions>& taskOptions)
{
    ocra::TaskConstructionManager factory(model, controller, taskOptions);
    return true;
}

void ControllerServer::setRegularizationTermWeights(double wDdq, double wTau, double wFc)
{
    if ( controllerType == WOCRA_CONTROLLER ) {
        std::dynamic_pointer_cast<wocra::WocraController>(controller)->setVariableMinimizationWeights(wDdq, wTau, wFc);
    }
}
