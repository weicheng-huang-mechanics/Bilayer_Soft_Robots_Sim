#ifndef DAMPINGFORCE_H
#define DAMPINGFORCE_H

#include "mechanics/base_force.h"

class dampingForce : public BaseForce
{
public:
	dampingForce(const std::shared_ptr<elasticPlate> m_plate, const std::shared_ptr<timeStepper> m_stepper, double m_viscosity);
	~dampingForce() override;

    void computeForceAndJacobian() override;
	void setFirstJacobian() override;

private:
	void computeFd();
	void computeJd();

    double viscosity;

    // support variables
    double ind;
    int index1;
    int index2;
    double edgeLength;
    double sectionArea;

};

#endif
