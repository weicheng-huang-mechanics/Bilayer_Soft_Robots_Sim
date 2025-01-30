#ifndef EXTERNALGRAVITYFORCE_H
#define EXTERNALGRAVITYFORCE_H

#include "mechanics/base_force.h"

class externalGravityForce : public BaseForce
{
public:
	externalGravityForce(const std::shared_ptr<elasticPlate> m_plate, const std::shared_ptr<timeStepper> m_stepper, Vector3d m_gVector);
	~externalGravityForce() override;

	void computeForceAndJacobian() override;
	void setFirstJacobian() override;

private:
	void computeFg();
	void computeJg();

	void setGravity();
	Vector3d gVector;	
    VectorXd massGravity;
};

#endif
