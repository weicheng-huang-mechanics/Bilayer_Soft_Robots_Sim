#include "elasticTwistingBound.h"
#include "elasticPlate.h"
#include "timeStepper.h"

elasticTwistingBound::elasticTwistingBound(const std::shared_ptr<elasticPlate> m_plate,
                                           const std::shared_ptr<timeStepper> m_stepper)
                                           : BaseForce(m_plate, m_stepper)
{	
    gradTwist = MatrixXd::Zero(plate->v_coupleBending.size(), 14);

    DDtwist.setZero(14,14);
    Jtt.setZero(14,14);
    gradTwistLocal.setZero(14);
    f.setZero(14);

    GJ = 1000 * plate->GJ;    

    localDOF = VectorXi::Zero(14);
}

elasticTwistingBound::~elasticTwistingBound() = default;

void elasticTwistingBound::computeForceAndJacobian()
{
	computeFt();
	computeJt();
}


void elasticTwistingBound::computeFt()
{
    for (int i = 0; i < plate->v_coupleBending.size(); i++)
    {
        deltam = plate->v_coupleBending[i].theta_2 - plate->v_coupleBending[i].theta_1;

        norm_e = plate->v_coupleBending[i].norm_1;
        norm_f = plate->v_coupleBending[i].norm_2;

        gradTwist.row(i).segment(0,3)  = -0.5 / norm_e * plate->v_coupleBending[i].kb;
        gradTwist.row(i).segment(11,3) = 0.5 / norm_f * plate->v_coupleBending[i].kb;
        gradTwist.row(i).segment(4,3)  = - gradTwist.row(i).segment(0,3);
        gradTwist.row(i).segment(7,3)  = - gradTwist.row(i).segment(11,3);

        gradTwist(i, 3)  = -1;
        gradTwist(i, 10) =  1;

        if ( plate->v_coupleBending[i].sign_1 < 0 )
        {
            gradTwist(i, 3) = - gradTwist(i, 3);
            
        }

        if ( plate->v_coupleBending[i].sign_2 < 0 )
        {
            gradTwist(i, 10) = - gradTwist(i, 10);
        }
        
    }

    for (int i = 0; i < plate->v_coupleBending.size(); i++)
    {
        deltam = plate->v_coupleBending[i].theta_2 - plate->v_coupleBending[i].theta_1;

        value = GJ / plate->v_coupleBending[i].voroniLength * (deltam + plate->v_coupleBending[i].refTwist);
        
        f = - value * gradTwist.row(i);
        localDOF = plate->v_coupleBending[i].arrayNum;

        for (int k = 0; k < 14; k++)
		{
			ind = localDOF(k);
			stepper->addForce(ind, -f[k]); // subtracting elastic force
		}
    }
}

void elasticTwistingBound::computeJt()
{
    for (int i = 0; i < plate->v_coupleBending.size(); i++)
    {
        norm_e = plate->v_coupleBending[i].norm_1;
        norm_f = plate->v_coupleBending[i].norm_2;
        te = plate->v_coupleBending[i].t_1;
        tf = plate->v_coupleBending[i].t_2;

        norm2_e=norm_e*norm_e;
        norm2_f=norm_f*norm_f;

        kbLocal = plate->v_coupleBending[i].kb;

        chi=1.0+te.dot(tf);
        tilde_t=(te+tf)/chi;

        crossMat(te,teMatrix);

        D2mDe2 = -0.25 / norm2_e * (kbLocal * (te+tilde_t).transpose()
            + (te+tilde_t) * kbLocal.transpose());
        D2mDf2 = -0.25 / norm2_f * (kbLocal * (tf+tilde_t).transpose()
            + (tf+tilde_t) * kbLocal.transpose());
        D2mDeDf = 0.5  / (norm_e*norm_f) * (2.0 / chi * teMatrix
            - kbLocal*tilde_t.transpose());
        D2mDfDe = D2mDeDf.transpose();


        //DDtwist.block(0,0,3,3) = D2mDe2;
        //DDtwist.block(0,4,3,3) =-D2mDe2 + D2mDeDf;
        //DDtwist.block(4,0,3,3) =-D2mDe2 + D2mDfDe;
        //DDtwist.block(4,4,3,3) = D2mDe2 - ( D2mDeDf + D2mDfDe ) + D2mDf2;
        //DDtwist.block(0,8,3,3) =-D2mDeDf;
        //DDtwist.block(8,0,3,3) =-D2mDfDe;
        //DDtwist.block(8,4,3,3) = D2mDfDe - D2mDf2;
        //DDtwist.block(4,8,3,3) = D2mDeDf - D2mDf2;
        //DDtwist.block(8,8,3,3) = D2mDf2;



        DDtwist.block(0,0,3,3)   = D2mDe2;
        DDtwist.block(4,4,3,3)   = D2mDe2;
        DDtwist.block(7,7,3,3)   = D2mDf2;
        DDtwist.block(11,11,3,3) = D2mDf2;

        DDtwist.block(0,4,3,3)   = - D2mDe2;
        DDtwist.block(4,0,3,3)   = - D2mDe2;
        DDtwist.block(7,11,3,3)  = - D2mDf2;
        DDtwist.block(11,7,3,3)  = - D2mDf2;

        DDtwist.block(0,7,3,3)   =   D2mDeDf;
        DDtwist.block(0,11,3,3)  = - D2mDeDf;
        DDtwist.block(7,0,3,3)   =   D2mDfDe;
        DDtwist.block(11,0,3,3)  = - D2mDfDe;

        DDtwist.block(4,7,3,3)   = - D2mDeDf;
        DDtwist.block(4,11,3,3)  =   D2mDeDf;
        DDtwist.block(7,4,3,3)   = - D2mDfDe;
        DDtwist.block(11,4,3,3)  =   D2mDfDe;



        gradTwistLocal = gradTwist.row(i);
        
        milen = -1 / plate->v_coupleBending[i].voroniLength;

        deltam = plate->v_coupleBending[i].theta_2 - plate->v_coupleBending[i].theta_1;

        Jtt = GJ * milen * ((deltam + plate->v_coupleBending[i].refTwist) 
            * DDtwist + gradTwistLocal * gradTwistLocal.transpose());

        localDOF = plate->v_coupleBending[i].arrayNum;

        for (int j=0;j<14;j++)
        {
            for (int k=0;k<14;k++)
            {
				ind1 = localDOF(j);
				ind2 = localDOF(k);
				stepper->addJacobian(ind1, ind2, - Jtt(k,j));
            }
        }
    }
}

// Utility
void elasticTwistingBound::crossMat(const Vector3d &a,Matrix3d &b)
{
	b<<0,-a(2),a(1),
	a(2),0,-a(0),
	-a(1),a(0),0;
}

void elasticTwistingBound::setFirstJacobian()
{
    for (int k = 0; k < plate->v_coupleBending.size(); k++)
    {
        localDOF = plate->v_coupleBending[k].arrayNum;

        for (int i = 0; i < 14; i++)
        {
            for (int j = 0; j < 14; j++)
            {
                stepper->addJacobian(localDOF(i), localDOF(j), 1);
            }
        }
    }
}