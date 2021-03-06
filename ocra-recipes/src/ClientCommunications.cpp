#include <ocra-recipes/ClientCommunications.h>


using namespace ocra_recipes;

int ClientCommunications::CONTROLLER_CLIENT_COUNT = 0;


ClientCommunications::ClientCommunications()
: inputCallback(*this)
{
    clientNumber = ++ClientCommunications::CONTROLLER_CLIENT_COUNT;

    while(yarp.exists(("/ControllerClient/"+ std::to_string(clientNumber) +"/rpc:o")))
    {
        ++clientNumber;
    }

    rpcClientPort_Name = "/ControllerClient/"+ std::to_string(clientNumber) +"/rpc:o";
    inputPort_Name = "/ControllerClient/"+ std::to_string(clientNumber) +":i";
}

ClientCommunications::~ClientCommunications()
{
    close();
    --ClientCommunications::CONTROLLER_CLIENT_COUNT;
}

bool ClientCommunications::open(double timeout, bool connectToTasks)
{
    rpcClientPort.open(rpcClientPort_Name.c_str());
    rpcClientPort.setReader(*this);
    inputPort.open(inputPort_Name.c_str());
    inputPort.setReader(inputCallback);

    bool isConOpen = openServerConnections(timeout);
    if(isConOpen && connectToTasks)
    {
        isConOpen &= openTaskConnections();

        if(isConOpen)
        {
            for(auto rpc_i : taskRpcClients)
            {
                yarp::os::Bottle message, reply;
                message.addInt(ocra::TASK_MESSAGE::GET_TYPE_AS_STRING);
                rpc_i.second->write(message, reply);
            }
        }else{
            yLog.error() << "Couldn't connect to the individual task ports.";
        }
    }
    return isConOpen;
}

std::vector<std::string> ClientCommunications::getTaskTypes()
{
    std::vector<std::string> retVec;
    for(auto rpc_i : taskRpcClients)
    {
        yarp::os::Bottle message, reply;
        message.addInt(ocra::TASK_MESSAGE::GET_TYPE_AS_STRING);
        rpc_i.second->write(message, reply);
        retVec.push_back(reply.get(0).asString());
    }
    return retVec;
}

bool ClientCommunications::close()
{
    rpcClientPort.close();
    inputPort.close();

    for(auto rpc_i : taskRpcClients)
    {
        rpc_i.second->close();
    }
    taskRpcClients.clear();
    return true;
}

void ClientCommunications::close(const std::string& taskName)
{
    if(taskRpcClients.find(taskName) != taskRpcClients.end())
    {
        taskRpcClients[taskName]->close();
        taskRpcClients.erase(taskName);
    }
}

bool ClientCommunications::read(yarp::os::ConnectionReader& connection)
{
    yarp::os::Bottle input;

    if (!input.read(connection)){
        return false;
    }
    else{
        parseMessage(input);
        return true;
    }
}

bool ClientCommunications::openServerConnections(double timeout)
{
    if (!yarp.checkNetwork()) {
        yLog.error() << "Yarp network isn't running.";
        return false;
    }
    else{
        bool connected = false;
        double timeDelayed = 0.0;
        double delayTime = 0.01;
        while(!connected && (timeDelayed < timeout))
        {
            connected = yarp.connect(rpcClientPort_Name.c_str(), "/ControllerServer/rpc:i");
            yarp::os::Time::delay(delayTime);
            timeDelayed += delayTime;
            if (timeDelayed>= timeout) {
                yLog.error() << "Could not connect to the ocra controller server. Are you sure it is running?";
            }
        }

        // connected = false;
        // timeDelayed = 0.0;
        // while(!connected && timeDelayed < timeout)
        // {
        //     connected = yarp.connect("/ControllerServer:o", inputPort_Name.c_str());
        //     yarp::os::Time::delay(delayTime);
        //     timeDelayed += delayTime;
        //     if (timeDelayed>= timeout) {
        //         yLog.error() << "Could not connect to the ocra controller port. Are you sure it is running?";
        //     }
        // }
        return connected;
    }
}
std::vector<std::string> ClientCommunications::getTaskPortNames()
{
    std::vector<std::string> portNameVec;
    yarp::os::Bottle message, reply;
    message.addInt(GET_TASK_PORT_LIST);
    rpcClientPort.write(message, reply);
    for(auto i=0; i<reply.size(); ++i)
    {
        portNameVec.push_back(reply.get(i).asString());
    }
    return portNameVec;
}

std::string ClientCommunications::getTaskPortName(const std::string& taskName)
{
    yarp::os::Bottle message, reply;
    message.addInt(GET_TASK_PORT_NAME);
    message.addString(taskName);
    rpcClientPort.write(message, reply);

    std::string portName = "";
    if (reply.size()>0) {
        portName = reply.get(0).asString();
    }
    return portName;
}

std::vector<std::string> ClientCommunications::getTaskNames()
{
    std::vector<std::string> nameVec;
    yarp::os::Bottle message, reply;
    message.addInt(GET_TASK_LIST);
    rpcClientPort.write(message, reply);
    for(auto i=0; i<reply.size(); ++i)
    {
        nameVec.push_back(reply.get(i).asString());
    }
    return nameVec;
}

std::shared_ptr<yarp::os::RpcClient> ClientCommunications::getTaskClient(const std::string& taskName)
{
    if(taskRpcClients.find(taskName) != taskRpcClients.end())
    {
        return taskRpcClients[taskName];
    }else{
        //TODO: return a proper null pointer.
    }
}


yarp::os::Bottle ClientCommunications::queryController(yarp::os::Bottle& requestBottle)
{
    yarp::os::Bottle reply;
    rpcClientPort.write(requestBottle, reply);
    return reply;
}

yarp::os::Bottle ClientCommunications::queryController(const SERVER_COMMUNICATIONS_MESSAGE request)
{
    yarp::os::Bottle requestBottle, reply;
    requestBottle.addInt(request);
    rpcClientPort.write(requestBottle, reply);
    return reply;
}

yarp::os::Bottle ClientCommunications::queryController(const std::vector<SERVER_COMMUNICATIONS_MESSAGE> requestVector)
{
    yarp::os::Bottle requestBottle, reply;
    for(auto request : requestVector){
        requestBottle.addInt(request);
    }
    rpcClientPort.write(requestBottle, reply);
    return reply;
}

bool ClientCommunications::openTaskConnections()
{
    std::vector<std::string> taskNames = getTaskNames();
    std::vector<std::string> taskPortNames = getTaskPortNames();

    bool taskConnected = taskNames.size() == taskPortNames.size();

    if(taskConnected)
    {
        for(auto i=0; i<taskPortNames.size(); ++i)
        {
            std::string tmpTaskPortName = "/ControllerClient/" + std::to_string(clientNumber) + "/" + taskNames[i] + ":o";
            taskRpcClients[taskNames[i]] = std::make_shared<yarp::os::RpcClient>();
            taskRpcClients[taskNames[i]]->open(tmpTaskPortName.c_str());
            taskConnected &= yarp.connect(tmpTaskPortName.c_str(), taskPortNames[i].c_str());
        }
    }else{
        yLog.error() << "The number of task ports and names does not match! Can't connect to task RPC ports.";
    }

    return taskConnected;
}

void ClientCommunications::parseMessage(yarp::os::Bottle& input)
{
    int btlSize = input.size();
    for (int i=0; i<btlSize;)
    {
        switch (input.get(i).asInt()) {
            case REMOVE_TASK_PORT:
                {
                    ++i;
                    std::cout << "Got message: REMOVE_TASK_PORT - " << input.get(i).asString() << std::endl;
                    close(input.get(i).asString());
                    ++i;
                }break;

            case HELP:
                {
                    ++i;
                    std::cout << "Got message: HELP." << std::endl;
                }break;

            default:
                {
                    ++i;
                    std::cout << "Got message: UNKNOWN." << std::endl;
                }break;

        }
    }
}


ClientCommunications::InputCallback::InputCallback(ClientCommunications& comsRef)
: coms(comsRef)
{
    //
}

bool ClientCommunications::InputCallback::read(yarp::os::ConnectionReader& connection)
{
    yarp::os::Bottle input;
    if (input.read(connection)){
        return coms.parseInput(input);
    }
    else{
        return false;
    }
}

bool ClientCommunications::parseInput(yarp::os::Bottle& input)
{
    int btlSize = input.size();
    for (auto i=0; i<btlSize; ++i) {
        switch (input.get(i).asInt()) {

            case REMOVE_TASK_PORT:
            {
                ++i;
                yLog.info() << "Removing task " << input.get(i).asString();
                close(input.get(i).asString());
            }break;

            default:
            break;
        }
    }
    return true;
}
