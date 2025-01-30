#ifndef CONSTRAINTFORCE_H
#define CONSTRAINTFORCE_H

#include "mechanics/base_force.h"
#include "symbolicEquations.h"

class constrainedForce : public BaseForce
{
    public:
        constrainedForce(const std::shared_ptr<elasticPlate> m_plate, const std::shared_ptr<timeStepper> m_stepper, double m_K, double h);
        ~constrainedForce() override;

        void computeForceAndJacobian() override;
        void setFirstJacobian() override;        

    private:

        shared_ptr<symbolicEquations> sym_eqs;

        Vector<double, 40> forces_input;

        Vector3d x1s, x1e, x2s, x2e;

        double theta_1, theta_2;

        Vector3d t1, t2, t10, t20;
        
        Vector3d x1s0, x1e0, x2s0, x2e0;
        Vector3d d1_10, d1_20, d2_10, d2_20;

        void prepareInputs();

        double K;
        double h;

        Vector<double, 14> forces;
        Matrix<double, 14, 14> jacobians;
};


#endif