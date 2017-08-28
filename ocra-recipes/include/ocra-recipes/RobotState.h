#ifndef ROBOT_STATE_H
#define ROBOT_STATE_H

#include <yarp/os/Portable.h>
#include <yarp/os/ConnectionWriter.h>
#include <yarp/os/ConnectionReader.h>
#include <yarp/os/PortReader.h>

#include <Eigen/Dense>
#include <Eigen/Lgsm>

#include <ocra/control/Model.h>
#include <ocra/util/Macros.h>
#include <ocra/util/ErrorsHelper.h>
#include <yarp/os/Bottle.h>


namespace ocra_recipes
{

/*! \class RobotState
 *  \brief A portable class for sending robot state information over yarp.
 *
 *  blah.
 */
class RobotState : public yarp::os::Portable
{
DEFINE_CLASS_POINTER_TYPEDEFS(RobotState)
public:
    RobotState ();
    RobotState (const int numberOfDoF);
    virtual ~RobotState ();

    virtual bool write(yarp::os::ConnectionWriter& connection);
    virtual bool read(yarp::os::ConnectionReader& connection);

    friend std::ostream& operator<<(std::ostream &out, const RobotState& state)
    {
//         std::cout << "q \t | \t qd" << std::endl;
//         for(auto i=0; i<state.q.size(); ++i)
//         {
//             out << state.q(i) << "\t | \t" << state.qd(i) << std::endl;
//         }
// 
//         out << "x" << " " << "y" << " " << "z" << " " << "qx" << " " << "qy" << " " << "qz" << " " << "qw" << std::endl;
//         out << state.H_root.x() << " ";
//         out << state.H_root.y() << " ";
//         out << state.H_root.z() << " ";
//         out << state.H_root.qx() << " ";
//         out << state.H_root.qy() << " ";
//         out << state.H_root.qz() << " ";
//         out << state.H_root.qw() << std::endl;
        out << "Homogeneous matrix LGSM: \n" << state.H_root.toHomogeneousMatrix() << std::endl;
        out << "Homogeneous matrix KDL \n" << state.H_rootKDL << std::endl;
        out << "rx" << " " << "ry" << " " << "rz" << " " << "rx" << " " << "ry" << " " << "rz" << std::endl;
        // out << state.T_root.rx() << " ";
        // out << state.T_root.ry() << " ";
        // out << state.T_root.rz() << " ";
        // out << state.T_root.vx() << " ";
        // out << state.T_root.vy() << " ";
        // out << state.T_root.vz() << std::endl;
        out << "Root twist LGSM \n" << state.T_root.transpose() << std::endl;
        out << "Root twist KDL \n" << state.T_rootKDL << std::endl;

        return out;
    }


public:
    Eigen::VectorXd              q;
    Eigen::VectorXd             qd;
    Eigen::Displacementd    H_root;
    Eigen::Twistd           T_root;
    // KDL Versions
    KDL::Frame              H_rootKDL;
    KDL::Twist              T_rootKDL;

private:
    int nDoF;
};



/*! \class StateListener
 *  \brief A callback for a port listing to robot states.
 *
 *  Blah.
 */
class StateListener : public yarp::os::PortReader {
private:
    std::shared_ptr<ocra::Model> model;

public:
    StateListener();
    StateListener (std::shared_ptr<ocra::Model> modelPtr);
    virtual ~StateListener ();

    bool read(yarp::os::ConnectionReader& connection);
};

} /* ocra_recipes */
#endif
