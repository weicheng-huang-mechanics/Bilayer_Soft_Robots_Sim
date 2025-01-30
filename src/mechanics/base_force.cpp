#include "base_force.h"

BaseForce::BaseForce(const std::shared_ptr<elasticPlate>& m_plate, const std::shared_ptr<timeStepper>& m_stepper)
    : plate(m_plate), stepper(m_stepper) {
}

BaseForce::~BaseForce() = default;

// void BaseForce::setTimeStepper(std::shared_ptr<BaseTimeStepper> m_stepper) {
//     stepper = m_stepper;
// }