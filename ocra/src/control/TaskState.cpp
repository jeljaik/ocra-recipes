#include <ocra/control/TaskState.h>

using namespace ocra;

TaskState::TaskState()
: containsPosition(false)
, containsVelocity(false)
, containsAcceleration(false)
, containsQ(false)
, containsQd(false)
, containsQdd(false)
, containsTorque(false)
, containsWrench(false)
, containsPositionKDL(false)
, containsVelocityKDL(false)
, containsAccelerationKDL(false)
, containsWrenchKDL(false)
{
    position = Eigen::Displacementd::Identity();
    velocity = Eigen::Twistd(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    acceleration = Eigen::Twistd(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    positionKDL = KDL::Frame();
    velocityKDL = KDL::Twist();
    accelerationKDL = KDL::Twist();
    q = Eigen::VectorXd::Zero(1);
    qd = Eigen::VectorXd::Zero(1);
    qdd = Eigen::VectorXd::Zero(1);
    torque = Eigen::VectorXd::Zero(1);

    wrench = Eigen::Wrenchd(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    wrenchKDL = KDL::Wrench();
}
TaskState::~TaskState()
{
    // Do nothing
}

Eigen::Displacementd TaskState::getPosition() const
{
    return this->position;
}

KDL::Frame TaskState::getPositionKDL() const
{
    return this->positionKDL;
}

Eigen::Twistd TaskState::getVelocity() const
{
    return this->velocity;
}

KDL::Twist TaskState::getVelocityKDL() const
{
    return this->velocityKDL;
}

Eigen::Twistd TaskState::getAcceleration() const
{
    return this->acceleration;
}

KDL::Twist TaskState::getAccelerationKDL() const
{
    return this->accelerationKDL;
}

Eigen::VectorXd TaskState::getQ() const
{
    return this->q;
}


Eigen::VectorXd TaskState::getQd() const
{
    return this->qd;
}


Eigen::VectorXd TaskState::getQdd() const
{
    return this->qdd;
}


Eigen::VectorXd TaskState::getTorque() const
{
    return this->torque;
}


Eigen::Wrenchd TaskState::getWrench() const
{
    return this->wrench;
}

KDL::Wrench TaskState::getWrenchKDL() const
{
    return this->wrenchKDL;
}


void TaskState::setPosition(const Eigen::Displacementd& newPosition)
{
    this->position = newPosition;
    this->containsPosition = true;
}

void TaskState::setPositionKDL(const KDL::Frame& newPosition)
{
    this->positionKDL = newPosition;
    this->containsPositionKDL = true;
}


void TaskState::setVelocity(const Eigen::Twistd& newVelocity)
{
    this->velocity = newVelocity;
    this->containsVelocity = true;
}

void TaskState::setVelocityKDL(const KDL::Twist& newVelocity)
{
    this->velocityKDL = newVelocity;
    this->containsVelocityKDL = true;
}


void TaskState::setAcceleration(const Eigen::Twistd& newAcceleration)
{
    this->acceleration = newAcceleration;
    this->containsAcceleration = true;
}

void TaskState::setAccelerationKDL(const KDL::Twist& newAcceleration)
{
    this->accelerationKDL = newAcceleration;
    this->containsAccelerationKDL = true;
}


void TaskState::setQ(const Eigen::VectorXd& newQ)
{
    this->q = newQ;
    this->containsQ = true;
}


void TaskState::setQd(const Eigen::VectorXd& newQd)
{
    this->qd = newQd;
    this->containsQd = true;
}


void TaskState::setQdd(const Eigen::VectorXd& newQdd)
{
    this->qdd = newQdd;
    this->containsQdd = true;
}


void TaskState::setTorque(const Eigen::VectorXd& newTorque)
{
    this->torque = newTorque;
    this->containsTorque = true;
}


void TaskState::setWrench(const Eigen::Wrenchd& newWrench)
{
    this->wrench = newWrench;
    this->containsWrench = true;
}

void TaskState::setWrenchKDL(const KDL::Wrench& newWrench)
{
    this->wrenchKDL = newWrench;
    this->containsWrenchKDL = true;
}

bool TaskState::hasPosition() const
{
    return this->containsPosition;
}
bool TaskState::hasPositionKDL() const
{
    return this->containsPositionKDL;
}

bool TaskState::hasVelocity() const
{
    return this->containsVelocity;
}
bool TaskState::hasVelocityKDL() const
{
    return this->containsVelocityKDL;
}

bool TaskState::hasAcceleration() const
{
    return this->containsAcceleration;
}
bool TaskState::hasAccelerationKDL() const
{
    return this->containsAccelerationKDL;
}

bool TaskState::hasQ() const
{
    return this->containsQ;
}

bool TaskState::hasQd() const
{
    return this->containsQd;
}

bool TaskState::hasQdd() const
{
    return this->containsQdd;
}

bool TaskState::hasTorque() const
{
    return this->containsTorque;
}

bool TaskState::hasWrench() const
{
    return this->containsWrench;
}
bool TaskState::hasWrenchKDL() const
{
    return this->containsWrenchKDL;
}

bool TaskState::extractFromBottle(const yarp::os::Bottle& bottle, int& sizeOfState)
{
    int i = 0;
    if (bottle.get(i).asInt() == TASK_STATE_BOTTLE)
    {
#ifdef OCRA_USES_KDL
        this->containsPositionKDL = bottle.get(++i).asBool();
        this->containsVelocityKDL = bottle.get(++i).asBool();
        this->containsAccelerationKDL = bottle.get(++i).asBool();
#else
        this->containsPosition = bottle.get(++i).asBool();
        this->containsVelocity = bottle.get(++i).asBool();
        this->containsAcceleration = bottle.get(++i).asBool();
#endif
        this->containsQ = bottle.get(++i).asBool();
        this->containsQd = bottle.get(++i).asBool();
        this->containsQdd = bottle.get(++i).asBool();
        this->containsTorque = bottle.get(++i).asBool();
#ifdef OCRA_USES_KDL
        this->containsWrenchKDL = bottle.get(++i).asBool();
#else
        this->containsWrench = bottle.get(++i).asBool();
#endif

        int indexesToSkip;

#ifdef OCRA_USES_KDL
        if (this->hasPositionKDL()) {
            this->setPositionKDL( util::pourBottleIntoFrame(util::trimBottle(bottle, i+1), indexesToSkip) );
            i += indexesToSkip;
        }
        if (this->hasVelocityKDL()) {
            this->setVelocityKDL( util::pourBottleIntoTwist(util::trimBottle(bottle, i+1), indexesToSkip) );
            i += indexesToSkip;
        }
        if (this->hasAccelerationKDL()) {
            this->setAccelerationKDL( util::pourBottleIntoTwist(util::trimBottle(bottle, i+1), indexesToSkip) );
            i += indexesToSkip;
        }
#else
        if (this->hasPosition()) {
            this->setPosition( util::pourBottleIntoDisplacementd(util::trimBottle(bottle, i+1), indexesToSkip) );
            i += indexesToSkip;
        }
        if (this->hasVelocity()) {
            this->setVelocity( util::pourBottleIntoTwistd(util::trimBottle(bottle, i+1), indexesToSkip) );
            i += indexesToSkip;
        }
        if (this->hasAcceleration()) {
            this->setAcceleration( util::pourBottleIntoTwistd(util::trimBottle(bottle, i+1), indexesToSkip) );
            i += indexesToSkip;
        }
#endif
        if (this->hasQ()) {
            this->setQ( util::pourBottleIntoEigenVector(util::trimBottle(bottle, i+1), indexesToSkip) );
            i += indexesToSkip;
        }
        if (this->hasQd()) {
            this->setQd( util::pourBottleIntoEigenVector(util::trimBottle(bottle, i+1), indexesToSkip) );
            i += indexesToSkip;
        }
        if (this->hasQdd()) {
            this->setQdd( util::pourBottleIntoEigenVector(util::trimBottle(bottle, i+1), indexesToSkip) );
            i += indexesToSkip;
        }
        if (this->hasTorque()) {
            this->setTorque( util::pourBottleIntoEigenVector(util::trimBottle(bottle, i+1), indexesToSkip) );
            i += indexesToSkip;
        }
#ifdef OCRA_USES_KDL
        if (this->hasWrenchKDL()) {
            this->setWrenchKDL( util::pourBottleIntoWrench(util::trimBottle(bottle, i+1), indexesToSkip) );
            i += indexesToSkip;
        }
#else
        if (this->hasWrench()) {
            this->setWrench( util::pourBottleIntoWrenchd(util::trimBottle(bottle, i+1), indexesToSkip) );
            i += indexesToSkip;
        }
#endif


        sizeOfState = i;
        return true;
    }
    return false;

}

void TaskState::putIntoBottle(yarp::os::Bottle& bottle) const
{
    bottle.addInt(TASK_STATE_BOTTLE);

#ifdef OCRA_USES_KDL
    bottle.addInt(this->hasPositionKDL());
    bottle.addInt(this->hasVelocityKDL());
    bottle.addInt(this->hasAccelerationKDL());
#else
    bottle.addInt(this->hasPosition());
    bottle.addInt(this->hasVelocity());
    bottle.addInt(this->hasAcceleration());
#endif
    
    bottle.addInt(this->hasQ());
    bottle.addInt(this->hasQd());
    bottle.addInt(this->hasQdd());
    bottle.addInt(this->hasTorque());
    
#ifdef OCRA_USES_KDL
    bottle.addInt(this->hasWrenchKDL());
#else
    bottle.addInt(this->hasWrench());
#endif

#ifdef OCRA_USES_KDL
    if (this->hasPositionKDL()) {
        util::pourFrameIntoBottle(this->getPositionKDL(), bottle);
    }
    if (this->hasVelocityKDL()) {
        util::pourTwistIntoBottle(this->getVelocityKDL(), bottle);
    }
    if (this->hasAccelerationKDL()) {
        util::pourTwistIntoBottle(this->getAccelerationKDL(), bottle);
    }
#else
    if (this->hasPosition()) {
        util::pourDisplacementdIntoBottle(this->getPosition(), bottle);
    }
    if (this->hasVelocity()) {
        util::pourTwistdIntoBottle(this->getVelocity(), bottle);
    }
    if (this->hasAcceleration()) {
        util::pourTwistdIntoBottle(this->getAcceleration(), bottle);
    }
#endif
    if (this->hasQ()) {
        util::pourEigenVectorIntoBottle(this->getQ(), bottle);
    }
    if (this->hasQd()) {
        util::pourEigenVectorIntoBottle(this->getQd(), bottle);
    }
    if (this->hasQdd()) {
        util::pourEigenVectorIntoBottle(this->getQdd(), bottle);
    }
    if (this->hasTorque()) {
        util::pourEigenVectorIntoBottle(this->getTorque(), bottle);
    }
#ifdef OCRA_USES_KDL
    if (this->hasWrenchKDL()) {
        util::pourWrenchIntoBottle(this->getWrenchKDL(), bottle);
    }
#else
    if (this->hasWrench()) {
        util::pourWrenchdIntoBottle(this->getWrench(), bottle);
    }
#endif


}
