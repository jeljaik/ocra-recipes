#include "ocra/control/Trajectory/TimeOptimalTrajectory.h"
#include <math.h>

namespace ocra
{

TimeOptimalTrajectory::~TimeOptimalTrajectory()
{
    delete gt_trajectory;
}


Eigen::MatrixXd TimeOptimalTrajectory::getDesiredValues(double _time)
{
    if (startTrigger)
    {
        startTrigger = false;
        t0 = _time;
    }

    Eigen::MatrixXd desiredValue = Eigen::MatrixXd::Zero(nDoF, TRAJ_DIM);

    double t = _time - t0;

    if (t <= duration) {
        desiredValue.col(POS_INDEX) = gt_trajectory->getPosition(t);
        desiredValue.col(VEL_INDEX) = gt_trajectory->getVelocity(t);
    } else {
        desiredValue.col(POS_INDEX) = gt_trajectory->getPosition(duration);
        desiredValue.col(VEL_INDEX) = gt_trajectory->getVelocity(duration);
    }
    return desiredValue;
}

void TimeOptimalTrajectory::initializeTrajectory()
{
    duration = 0.0;
    /* For reaching */
    // maximumVelocityVector = Eigen::VectorXd::Constant(nDoF, 0.05);
    // maximumAccelerationVector = Eigen::VectorXd::Constant(nDoF, 0.05);

    /* For reaching on real robot */
    // maximumVelocityVector = Eigen::VectorXd::Constant(nDoF, 0.05); // will work for original case
    // maximumAccelerationVector = Eigen::VectorXd::Constant(nDoF, 0.05); // will work for original case
    // maximumVelocityVector = Eigen::VectorXd::Constant(nDoF, 0.01);
	// maximumAccelerationVector = Eigen::VectorXd::Constant(nDoF, 0.01);


    /* For standing up */
    maximumVelocityVector = Eigen::VectorXd::Constant(nDoF, 0.1);
    maximumAccelerationVector = Eigen::VectorXd::Constant(nDoF, 0.1);
    maxDeviation = 0.1;
    timeStep = 0.01;

    recalculateTrajectory();
}

double TimeOptimalTrajectory::getDuration()
{
    return duration;
}

void TimeOptimalTrajectory::recalculateTrajectory()
{
    gttraj::Path path = gttraj::Path(waypointList, maxDeviation);
    gt_trajectory = new gttraj::Trajectory(path, maximumVelocityVector, maximumAccelerationVector, timeStep);
    if(gt_trajectory->isValid()) {
		duration = gt_trajectory->getDuration();
    }
	else {
		OCRA_WARNING("Trajectory generation failed.")
	}
}



} //namespace ocra
