#include <ocra-recipes/RobotState.h>

using namespace ocra_recipes;


RobotState::RobotState()
{
}

RobotState::RobotState(const int numberOfDoF)
: nDoF(numberOfDoF)
{
    q.resize(nDoF);
    qd.resize(nDoF);
}

RobotState::~RobotState()
{
}

bool RobotState::write(yarp::os::ConnectionWriter& connection)
{
    // Two lists will be transmitted.
    connection.appendInt(BOTTLE_TAG_LIST);
    connection.appendInt(2);
    // The first list will contain just an integer (DOF))
    connection.appendInt(BOTTLE_TAG_INT);
    connection.appendInt(q.size());
    // The second list will contain the state of the robot
    connection.appendInt(BOTTLE_TAG_LIST + BOTTLE_TAG_DOUBLE);
    connection.appendInt(q.size() + qd.size() + 13);
    for(auto i=0; i<q.size(); ++i)
    {
        connection.appendDouble(q(i));
        connection.appendDouble(qd(i));
    }
    connection.appendDouble(H_root.x());
    connection.appendDouble(H_root.y());
    connection.appendDouble(H_root.z());
    connection.appendDouble(H_root.qx());
    connection.appendDouble(H_root.qy());
    connection.appendDouble(H_root.qz());
    connection.appendDouble(H_root.qw());

    connection.appendDouble(T_root.rx());
    connection.appendDouble(T_root.ry());
    connection.appendDouble(T_root.rz());
    connection.appendDouble(T_root.vx());
    connection.appendDouble(T_root.vy());
    connection.appendDouble(T_root.vz());
    
    // If someone connects in text mode, show something readable
    connection.convertTextMode();
    
    return !connection.isError();
}

bool RobotState::read(yarp::os::ConnectionReader& connection)
{
    // Auto-convert text mode interaction
    connection.convertTextMode();
    
    if ( connection.expectInt() != BOTTLE_TAG_LIST || connection.expectInt() != 2) {
        OCRA_ERROR("Received malformed data. Expected two lists with the following structure: (DOF)(ROBOT STATE)");
        return false;
    }
    
    // Reading DOF
    if ( connection.expectInt() != BOTTLE_TAG_INT ) {
        OCRA_ERROR("Received malformed data. Expected one integer for DOF)");
        return false;
    }
    this->nDoF = connection.expectInt();
    this->q.resize(nDoF);
    this->qd.resize(nDoF);

    if ( connection.expectInt() != BOTTLE_TAG_LIST + BOTTLE_TAG_DOUBLE || connection.expectInt()!= q.size() + qd.size() + 13 ) {
        OCRA_ERROR("Received a list with less data than expected");
        return false;
    }

    for(auto i=0; i<this->nDoF; ++i)
    {
        this->q(i) = connection.expectDouble();
        this->qd(i) = connection.expectDouble();
    }
    this->H_root.x() = connection.expectDouble();
    this->H_root.y() = connection.expectDouble();
    this->H_root.z() = connection.expectDouble();
    this->H_root.qx() = connection.expectDouble();
    this->H_root.qy() = connection.expectDouble();
    this->H_root.qz() = connection.expectDouble();
    this->H_root.qw() = connection.expectDouble();

    this->T_root.rx() = connection.expectDouble();
    this->T_root.ry() = connection.expectDouble();
    this->T_root.rz() = connection.expectDouble();
    this->T_root.vx() = connection.expectDouble();
    this->T_root.vy() = connection.expectDouble();
    this->T_root.vz() = connection.expectDouble();
    
    return !connection.isError();
}


/**************************************************************************************************
                                    StateListener Class
**************************************************************************************************/
StateListener::StateListener()
{
}

StateListener::~StateListener()
{
}

StateListener::StateListener(std::shared_ptr<ocra::Model> modelPtr)
: model(modelPtr)
{
    //do nothing
}

bool StateListener::read(yarp::os::ConnectionReader& connection)
{

    RobotState state;

    if (!state.read(connection)){
        OCRA_ERROR("Couldn't read state: " << state << std::endl);
        return false;
    }
    else{
        model->setState(state.H_root, state.q, state.T_root, state.qd);
        return true;
    }
}
