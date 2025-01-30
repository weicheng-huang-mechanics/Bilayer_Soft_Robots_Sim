#ifndef BASE_FORCE_H
#define BASE_FORCE_H

#include "eigenIncludes.h"

class timeStepper;
class elasticPlate;

class BaseForce
{
  public:
    explicit BaseForce(const std::shared_ptr<elasticPlate>& m_plate, const std::shared_ptr<timeStepper>& m_stepper);
    virtual ~BaseForce() = 0;

    virtual void computeForceAndJacobian() = 0;
    virtual	void setFirstJacobian() = 0;
    // void setTimeStepper(std::shared_ptr<BaseTimeStepper> m_stepper);

  protected:
    std::shared_ptr<elasticPlate> plate = nullptr;
    std::shared_ptr<timeStepper> stepper = nullptr;
};

#endif  // BASE_FORCE_H