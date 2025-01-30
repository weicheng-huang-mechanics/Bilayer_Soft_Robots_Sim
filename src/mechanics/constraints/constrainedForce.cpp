#include "constrainedForce.h"
#include "elasticPlate.h"
#include "timeStepper.h"


constrainedForce::constrainedForce(const std::shared_ptr<elasticPlate> m_plate, 
                                   const std::shared_ptr<timeStepper> m_stepper, double m_K, double m_h)
                                   : BaseForce(m_plate, m_stepper), K(m_K), h(m_h)
{
    sym_eqs = make_shared<symbolicEquations>();
    sym_eqs->generatePenaltyFunctions();

    forces_input(38) = K;
    forces_input(39) = h;
}

constrainedForce::~constrainedForce() = default;


void constrainedForce::computeForceAndJacobian()
{
    for (int k = 0; k < plate->coupleBend.size(); k++)
    {

        int idx1 = plate->coupleBend[k](0);
        int idx2 = plate->coupleBend[k](1);

        int nv_1 = plate->v_edgeElement[idx1].nv_1;
        int nv_2 = plate->v_edgeElement[idx1].nv_2;

        x1s = plate->v_edgeElement[idx1].x_1;
        x1e = plate->v_edgeElement[idx1].x_2;
        theta_1 = plate->v_edgeElement[idx1].theta;

        x1s0 = plate->v_edgeElement[idx1].x_1_old;
        x1e0 = plate->v_edgeElement[idx1].x_2_old;

        d1_10 = plate->v_edgeElement[idx1].d1_old;
        d1_20 = plate->v_edgeElement[idx1].d2_old;

        int nv_11 = plate->v_edgeElement[idx2].nv_1;
        int nv_22 = plate->v_edgeElement[k + 99].nv_2;

        x2s = plate->v_edgeElement[k + 99].x_1;
        x2e = plate->v_edgeElement[k + 99].x_2;
        theta_2 = plate->v_edgeElement[k + 99].theta;

        x2s0 = plate->v_edgeElement[k + 99].x_1_old;
        x2e0 = plate->v_edgeElement[k + 99].x_2_old;

        d2_10 = plate->v_edgeElement[k + 99].d1_old;
        d2_20 = plate->v_edgeElement[k + 99].d2_old;

        t1 = (x1e - x1s) / (x1e - x1s).norm();
        t10 = (x1e0 - x1s0) / (x1e0 - x1s0).norm();
        t2 = (x2e - x2s) / (x2e - x2s).norm();
        t20 = (x2e0 - x2s0) / (x2e0 - x2s0).norm();

        prepareInputs();

        if (t1.cross(t10).norm() < 1.0e-12 && t2.cross(t20).norm() < 1.0e-12)
        {
            sym_eqs->E_grad_dx_func_3.call(forces.data(), forces_input.data());
            sym_eqs->E_hess_dx_func_3.call(jacobians.data(), forces_input.data());

        }
        else if (t1.cross(t10).norm() < 1.0e-12)
        {
            sym_eqs->E_grad_dx_func_1.call(forces.data(), forces_input.data());
            sym_eqs->E_hess_dx_func_1.call(jacobians.data(), forces_input.data()); 
        }
        else if (t2.cross(t20).norm() < 1.0e-12){
            sym_eqs->E_grad_dx_func_2.call(forces.data(), forces_input.data());
            sym_eqs->E_hess_dx_func_2.call(jacobians.data(), forces_input.data()); 

        }
        else{
            sym_eqs->E_grad_dx_func.call(forces.data(), forces_input.data());
            sym_eqs->E_hess_dx_func.call(jacobians.data(), forces_input.data()); 
        } 
    
        VectorXi localDOF(14);
        localDOF << 3 * nv_1, 3 * nv_1 + 1, 3 * nv_1 + 2, 3 * plate->nv + nv_1,
                    3 * nv_2, 3 * nv_2 + 1, 3 * nv_2 + 2,
                    3 * nv_11, 3 * nv_11 + 1, 3 * nv_11 + 2, 3 * plate->nv + nv_11,
                    3 * nv_22, 3 * nv_22 + 1, 3 * nv_22 + 2;

        for(int i = 0; i < 14; i++)
		{
			int ind1 = localDOF(i);
			stepper->addForce(ind1, forces[i]);
		}

        for(int i = 0; i < 14; i++)
		{

			int ind1 = localDOF(i);         
            for (int j = 0; j < 14; j++)
            {
                int ind2 = localDOF(j);
			    stepper->addJacobian(ind1, ind2, jacobians(i, j));
            }
		}

    }
}



void constrainedForce::prepareInputs()
{
    forces_input(seq(0, 2)) = x1s;
    forces_input(3) = theta_1;
    forces_input(seq(4, 6)) = x1e;
    forces_input(seq(7, 9)) = x2s;
    forces_input(10) = theta_2;
    forces_input(seq(11, 13)) = x2e;
    forces_input(seq(14, 16)) = x1s0;
    forces_input(seq(17, 19)) = x1e0;
    forces_input(seq(20, 22)) = x2s0;
    forces_input(seq(23, 25)) = x2e0;
    forces_input(seq(26, 28)) = d1_10;
    forces_input(seq(29, 31)) = d1_20;
    forces_input(seq(32, 34)) = d2_10;
    forces_input(seq(35, 37)) = d2_20;
}

void constrainedForce::setFirstJacobian()
{
    ;
}