#include <ocra/control/TaskBuilders/CartesianTaskBuilder.h>

using namespace ocra;

CartesianTaskBuilder::CartesianTaskBuilder(const TaskBuilderOptions& taskOptions, Model::Ptr modelPtr)
: TaskBuilder(taskOptions, modelPtr)
{
    this->nDoF = utils::computeDimensionFor(this->options.axes);
}

CartesianTaskBuilder::~CartesianTaskBuilder()
{

}

Feature::Ptr CartesianTaskBuilder::buildFeature()
{
    std::string featFrameName = this->options.taskName + ".SegmentFrame";
    std::string featName = this->options.taskName + ".PositionFeature";
    std::string segmentName = this->model->SegmentName(this->options.segment);

    ControlFrame::Ptr featFrame =  std::make_shared<SegmentFrame>(featFrameName, *this->model, segmentName, this->options.offset);

    return std::make_shared<PositionFeature>(featName, featFrame, this->options.axes);
}

Feature::Ptr CartesianTaskBuilder::buildFeatureDesired()
{
    std::string featDesFrameName = this->options.taskName + ".TargetFrame";
    std::string featDesName = this->options.taskName + ".PositionFeatureDesired";

    ControlFrame::Ptr featDesFrame = std::make_shared<TargetFrame>(featDesFrameName, *this->model);

    return std::make_shared<PositionFeature>(featDesName, featDesFrame, this->options.axes);
}

void CartesianTaskBuilder::setTaskState()
{
    TaskState state;

    // Make sure the desired vector is the same size as the number of axes being controlled.
    if(this->options.desired.size() > 0){
        // Set the pose
#ifndef OCRA_USES_KDL
        state.setPosition(util::eigenVectorToDisplacementd(this->options.desired));
#else
        KDL::Frame tmpFrame;
        util::eigenVectorToKDLFrame(this->options.desired, tmpFrame);
        state.setPositionKDL(tmpFrame);
#endif 
    }else{
        // If the desired position was not given then just get it from the current state of the task.
#ifndef OCRA_USES_KDL
        state.setPosition(this->task->getTaskState().getPosition());
#else 
        state.setPositionKDL(this->task->getTaskState().getPositionKDL());
#endif
    }

    // Set velocity and acceleration to zero.
#ifndef OCRA_USES_KDL
    state.setVelocity(Eigen::Twistd::Zero());
    state.setAcceleration(Eigen::Twistd::Zero());
#else
    state.setVelocityKDL(KDL::Twist::Zero());
    state.setAccelerationKDL(KDL::Twist::Zero());
#endif
    
    this->task->setDesiredTaskStateDirect(state);
}

void CartesianTaskBuilder::setTaskType()
{
    this->task->setTaskType(Task::ACCELERATIONTASK);
    this->task->setMetaTaskType(Task::POSITION);
}
