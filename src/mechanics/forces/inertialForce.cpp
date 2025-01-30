#include "inertialForce.h"
#include "elasticPlate.h"
#include "timeStepper.h"

inertialForce::inertialForce(const std::shared_ptr<elasticPlate> m_plate, 
							 const std::shared_ptr<timeStepper> m_stepper)
							 : BaseForce(m_plate, m_stepper)
{
}

inertialForce::~inertialForce() = default;

void inertialForce::computeForceAndJacobian()
{
	computeFi();
	computeJi();
}

void inertialForce::computeFi()
{
	for (int i = 0; i < plate->ndof; i++)
	{
		f = plate->massArray[i] * (plate->x[i] - plate->x0[i]) / ((plate->dt) *(plate->dt))
				- (plate->massArray[i] * plate->u[i])/(plate->dt);

		stepper->addForce(i, f);
	}
}

void inertialForce::computeJi()
{
	for (int i=0; i<plate->ndof; i++)
    {
		jac = plate->massArray(i)/ ((plate->dt) *(plate->dt));
		stepper->addJacobian(i, i, jac);
	}
}

void inertialForce::setFirstJacobian()
{
	for (int i=0; i<plate->ndof; i++)
	{
		stepper->addJacobian(i, i, 1);
	}
}
