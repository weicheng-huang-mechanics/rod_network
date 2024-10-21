#include "elasticBendingForce.h"

elasticBendingForce::elasticBendingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
    stepper = &m_stepper;

	Id3<<1,0,0,
         0,1,0,
         0,0,1;

	double EI = plate->EI;
    EIMat<<EI,0,
           0,EI;

	int bendingNum = plate->bendingNum;
    gradKappa1 = MatrixXd::Zero(bendingNum,11);
    gradKappa2 = MatrixXd::Zero(bendingNum,11);
    relevantPart = MatrixXd::Zero(11, 2);;
    DDkappa1 = MatrixXd::Zero(11,11);
    DDkappa2 = MatrixXd::Zero(11,11);
    Jbb = MatrixXd::Zero(11,11);
    
    D2kappa1De2.setZero(3,3);
    D2kappa1Df2.setZero(3,3);
    D2kappa1DeDf.setZero(3,3);
    D2kappa2De2.setZero(3,3);
    D2kappa2Df2.setZero(3,3);
    D2kappa2DeDf.setZero(3,3);
    kappa11 = VectorXd::Zero(bendingNum);
    kappa22 = VectorXd::Zero(bendingNum);
    f = VectorXd::Zero(11);

    localDOF = VectorXi::Zero(11);
}

elasticBendingForce::~elasticBendingForce()
{
	;
}

void elasticBendingForce::computeFb()
{
    
    for (int i = 0; i < plate->bendingNum; i++)
    {
        norm_e = plate->v_bendingElement[i].norm_1;
        norm_f = plate->v_bendingElement[i].norm_2;
        te = plate->v_bendingElement[i].t_1;
        tf = plate->v_bendingElement[i].t_2;
        d1e = plate->v_bendingElement[i].m_11;
        d2e = plate->v_bendingElement[i].m_12;
        d1f = plate->v_bendingElement[i].m_21;
        d2f = plate->v_bendingElement[i].m_22;

        chi = 1.0 + te.dot(tf);
        tilde_t = (te+tf)/chi;
        tilde_d1 = (d1e+d1f)/chi;
        tilde_d2 = (d2e+d2f)/chi;

        kappa1 = plate->v_bendingElement[i].kappa(0);
        kappa2 = plate->v_bendingElement[i].kappa(1);

        Dkappa1De = (1.0/norm_e)*(-kappa1*tilde_t + tf.cross(tilde_d2));
        Dkappa1Df = (1.0/norm_f)*(-kappa1*tilde_t - te.cross(tilde_d2));
        Dkappa2De = (1.0/norm_e)*(-kappa2*tilde_t - tf.cross(tilde_d1));
        Dkappa2Df = (1.0/norm_f)*(-kappa2*tilde_t + te.cross(tilde_d1));

        gradKappa1.row(i).segment(0,3)=-Dkappa1De;
        gradKappa1.row(i).segment(4,3)= Dkappa1De - Dkappa1Df;
        gradKappa1.row(i).segment(8,3)= Dkappa1Df;

        gradKappa2.row(i).segment(0,3)=-Dkappa2De;
        gradKappa2.row(i).segment(4,3)= Dkappa2De - Dkappa2Df;
        gradKappa2.row(i).segment(8,3)= Dkappa2Df;

        kbLocal = plate->v_bendingElement[i].kb;

        gradKappa1(i,3)=-0.5*kbLocal.dot(d1e);
        gradKappa1(i,7)=-0.5*kbLocal.dot(d1f);
        gradKappa2(i,3)=-0.5*kbLocal.dot(d2e);
        gradKappa2(i,7)=-0.5*kbLocal.dot(d2f);

        if ( plate->v_bendingElement[i].sign_1 < 0 )
        {
            gradKappa1(i,3) = - gradKappa1(i,3);
            gradKappa2(i,3) = - gradKappa2(i,3);
        }

        if ( plate->v_bendingElement[i].sign_2 < 0 )
        {
            gradKappa1(i,7) = - gradKappa1(i,7);
            gradKappa2(i,7) = - gradKappa2(i,7);
        }
    }

    for (int i = 0; i < plate->bendingNum; i++)
    {
        localDOF = plate->v_bendingElement[i].arrayNum;

        relevantPart.col(0) = gradKappa1.row(i);
        relevantPart.col(1) = gradKappa2.row(i);
        kappaL = plate->v_bendingElement[i].kappa - plate->v_bendingElement[i].kappaBar;
        f = - relevantPart * EIMat * kappaL / plate->v_bendingElement[i].voroniLength;

        for (int k = 0; k < 11; k++)
		{
			int ind = localDOF(k);
			stepper->addForce(ind, -f[k]); // subtracting elastic force
		}
    }
    
}

void elasticBendingForce::computeJb()
{
    
	for (int i = 0; i < plate->bendingNum; i++)
    {
        norm_e = plate->v_bendingElement[i].norm_1;
        norm_f = plate->v_bendingElement[i].norm_2;
        te = plate->v_bendingElement[i].t_1;
        tf = plate->v_bendingElement[i].t_2;
        d1e = plate->v_bendingElement[i].m_11;
        d2e = plate->v_bendingElement[i].m_12;
        d1f = plate->v_bendingElement[i].m_21;
        d2f = plate->v_bendingElement[i].m_22;

        norm2_e=norm_e*norm_e;
        norm2_f=norm_f*norm_f;

        chi=1.0+te.dot(tf);
        tilde_t=(te+tf)/chi;
        tilde_d1=(d1e+d1f)/chi;
        tilde_d2=(d2e+d2f)/chi;

        kappa1 = plate->v_bendingElement[i].kappa(0);
        kappa2 = plate->v_bendingElement[i].kappa(1);

        kbLocal = plate->v_bendingElement[i].kb;
		
        tt_o_tt = tilde_t * tilde_t.transpose();

        crossMat(tilde_d1,tilde_d1_3d);
        crossMat(tilde_d2,tilde_d2_3d);

        tmp = tf.cross(tilde_d2);
        tf_c_d2t_o_tt=tmp*tilde_t.transpose();
        tt_o_tf_c_d2t=tf_c_d2t_o_tt.transpose();
        kb_o_d2e=kbLocal*d2e.transpose();
        d2e_o_kb=kb_o_d2e.transpose();

        D2kappa1De2=
        1.0/norm2_e * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt - tt_o_tf_c_d2t)
        - kappa1 / (chi * norm2_e) * (Id3 - te * te.transpose())
        + 1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb);

        tmp=te.cross(tilde_d2);
        te_c_d2t_o_tt = tmp * tilde_t.transpose();
        tt_o_te_c_d2t= te_c_d2t_o_tt.transpose();
        kb_o_d2f= kbLocal * d2f.transpose();
        d2f_o_kb= kb_o_d2f.transpose();

        D2kappa1Df2=
        1.0 / norm2_f * (2* kappa1 * tt_o_tt + te_c_d2t_o_tt + tt_o_te_c_d2t)
        - kappa1 / (chi * norm2_f) * (Id3 - tf*tf.transpose())
        + 1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb);

        D2kappa1DeDf=
        - kappa1/(chi* norm_e * norm_f) *(Id3 + te*tf.transpose())
        + 1.0 / (norm_e*norm_f) * (2*kappa1 *tt_o_tt - tf_c_d2t_o_tt +
            tt_o_te_c_d2t - tilde_d2_3d);
        D2kappa1DfDe = D2kappa1DeDf.transpose();

        tmp = tf.cross(tilde_d1);
        tf_c_d1t_o_tt = tmp*tilde_t.transpose();
        tt_o_tf_c_d1t = tf_c_d1t_o_tt.transpose();
        kb_o_d1e = kbLocal * d1e.transpose();
        d1e_o_kb = kb_o_d1e.transpose();

        D2kappa2De2
        = 1.0 / norm2_e * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt + tt_o_tf_c_d1t)
        - kappa2 / (chi * norm2_e) * (Id3 - te*te.transpose())
        - 1.0 / (4.0 * norm2_e) * (kb_o_d1e + d1e_o_kb);

        tmp = te.cross(tilde_d1);
        Matrix3d te_c_d1t_o_tt = tmp*tilde_t.transpose();
        Matrix3d tt_o_te_c_d1t = te_c_d1t_o_tt.transpose();
        Matrix3d kb_o_d1f = kbLocal*d1f.transpose();
        Matrix3d d1f_o_kb = kb_o_d1f.transpose();

        D2kappa2Df2
        = 1.0 / norm2_f * (2 * kappa2 * tt_o_tt - te_c_d1t_o_tt - tt_o_te_c_d1t)
        - kappa2 / (chi * norm2_f) * (Id3 - tf*tf.transpose())
        - 1.0 / (4.0 * norm2_f) * (kb_o_d1f + d1f_o_kb);

        D2kappa2DeDf
        = -kappa2/(chi * norm_e * norm_f) * (Id3 + te*tf.transpose())
        + 1.0 / (norm_e*norm_f) * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt 
        - tt_o_te_c_d1t + tilde_d1_3d);

        D2kappa2DfDe = D2kappa2DeDf.transpose();

        D2kappa1Dthetae2 = -0.5 * kbLocal.dot(d2e);
        D2kappa1Dthetaf2 = -0.5 * kbLocal.dot(d2f);
        D2kappa2Dthetae2 =  0.5 * kbLocal.dot(d1e);
        D2kappa2Dthetaf2 =  0.5 * kbLocal.dot(d1f);

        D2kappa1DeDthetae
        = 1.0 / norm_e * ((0.5 * kbLocal.dot(d1e)) * tilde_t 
        - 1.0 / chi * (tf.cross(d1e)));
        D2kappa1DeDthetaf
        = 1.0 / norm_e * ((0.5 * kbLocal.dot(d1f)) * tilde_t
        - 1.0 / chi * (tf.cross(d1f)));
        D2kappa1DfDthetae
        = 1.0 / norm_f * ((0.5 * kbLocal.dot(d1e)) * tilde_t 
        + 1.0 / chi * (te.cross(d1e)));
        D2kappa1DfDthetaf
        = 1.0 / norm_f * ((0.5 * kbLocal.dot(d1f)) * tilde_t 
        + 1.0 / chi * (te.cross(d1f)));
        D2kappa2DeDthetae
        = 1.0 / norm_e * ((0.5 * kbLocal.dot(d2e)) * tilde_t
        - 1.0 / chi * (tf.cross(d2e)));
        D2kappa2DeDthetaf
        = 1.0 / norm_e * ((0.5 * kbLocal.dot(d2f)) * tilde_t 
        - 1.0 / chi * (tf.cross(d2f)));
        D2kappa2DfDthetae
        = 1.0 / norm_f * ((0.5 * kbLocal.dot(d2e)) * tilde_t
        + 1.0 / chi * (te.cross(d2e)));
        D2kappa2DfDthetaf
        = 1.0 / norm_f * ((0.5 * kbLocal.dot(d2f)) * tilde_t
        + 1.0 / chi * (te.cross(d2f)));

        DDkappa1.block(0,0,3,3) =   D2kappa1De2;
        DDkappa1.block(0,4,3,3) = - D2kappa1De2 + D2kappa1DeDf;
        DDkappa1.block(0,8,3,3) =               - D2kappa1DeDf;
        DDkappa1.block(4,0,3,3) = - D2kappa1De2                + D2kappa1DfDe;
        DDkappa1.block(4,4,3,3) =   D2kappa1De2 - D2kappa1DeDf - D2kappa1DfDe + D2kappa1Df2;
        DDkappa1.block(4,8,3,3) =                 D2kappa1DeDf                - D2kappa1Df2;
        DDkappa1.block(8,0,3,3) =                              - D2kappa1DfDe;
        DDkappa1.block(8,4,3,3) =                                D2kappa1DfDe - D2kappa1Df2;
        DDkappa1.block(8,8,3,3) =                                               D2kappa1Df2;

        DDkappa1(3, 3) =   D2kappa1Dthetae2;
        DDkappa1(7, 7) =   D2kappa1Dthetaf2;

        DDkappa1.col(3).segment(0,3) = - D2kappa1DeDthetae;
        DDkappa1.col(3).segment(4,3) =   D2kappa1DeDthetae - D2kappa1DfDthetae;
        DDkappa1.col(3).segment(8,3) =                       D2kappa1DfDthetae;
        DDkappa1.row(3).segment(0,3) =   DDkappa1.col(3).segment(0,3).transpose();
        DDkappa1.row(3).segment(4,3) =   DDkappa1.col(3).segment(4,3).transpose();
        DDkappa1.row(3).segment(8,3) =   DDkappa1.col(3).segment(8,3).transpose();

        DDkappa1.col(7).segment(0,3) = - D2kappa1DeDthetaf;
        DDkappa1.col(7).segment(4,3) =   D2kappa1DeDthetaf - D2kappa1DfDthetaf;
        DDkappa1.col(7).segment(8,3) =                       D2kappa1DfDthetaf;
        DDkappa1.row(7).segment(0,3) =   DDkappa1.col(7).segment(0,3).transpose();
        DDkappa1.row(7).segment(4,3) =   DDkappa1.col(7).segment(4,3).transpose();
        DDkappa1.row(7).segment(8,3) =   DDkappa1.col(7).segment(8,3).transpose();

        DDkappa2.block(0,0,3,3) =   D2kappa2De2;
        DDkappa2.block(0,4,3,3) = - D2kappa2De2 + D2kappa2DeDf;
        DDkappa2.block(0,8,3,3) =               - D2kappa2DeDf;
        DDkappa2.block(4,0,3,3) = - D2kappa2De2                + D2kappa2DfDe;
        DDkappa2.block(4,4,3,3) =   D2kappa2De2 - D2kappa2DeDf - D2kappa2DfDe + D2kappa2Df2;
        DDkappa2.block(4,8,3,3) =                 D2kappa2DeDf                - D2kappa2Df2;
        DDkappa2.block(8,0,3,3) =                              - D2kappa2DfDe;
        DDkappa2.block(8,4,3,3) =                                D2kappa2DfDe - D2kappa2Df2;
        DDkappa2.block(8,8,3,3) =                                               D2kappa2Df2;

        DDkappa2(3, 3)     = D2kappa2Dthetae2;
        DDkappa2(7, 7)     = D2kappa2Dthetaf2;

        DDkappa2.col(3).segment(0,3) = - D2kappa2DeDthetae;
        DDkappa2.col(3).segment(4,3) =   D2kappa2DeDthetae - D2kappa2DfDthetae;
        DDkappa2.col(3).segment(8,3) =                       D2kappa2DfDthetae;
        DDkappa2.row(3).segment(0,3) =   DDkappa2.col(3).segment(0,3).transpose();
        DDkappa2.row(3).segment(4,3) =   DDkappa2.col(3).segment(4,3).transpose();
        DDkappa2.row(3).segment(8,3) =   DDkappa2.col(3).segment(8,3).transpose();

        DDkappa2.col(7).segment(0,3) = - D2kappa2DeDthetaf;
        DDkappa2.col(7).segment(4,3) =   D2kappa2DeDthetaf - D2kappa2DfDthetaf;
        DDkappa2.col(7).segment(8,3) =                       D2kappa2DfDthetaf;
        DDkappa2.row(7).segment(0,3) =   DDkappa2.col(7).segment(0,3).transpose();
        DDkappa2.row(7).segment(4,3) =   DDkappa2.col(7).segment(4,3).transpose();
        DDkappa2.row(7).segment(8,3) =   DDkappa2.col(7).segment(8,3).transpose();

        if ( plate->v_bendingElement[i].sign_1 < 0)
        {
            DDkappa1.col(3).segment(0,3) = - DDkappa1.col(3).segment(0,3);
            DDkappa1.col(3).segment(4,3) = - DDkappa1.col(3).segment(4,3);
            DDkappa1.col(3).segment(8,3) = - DDkappa1.col(3).segment(8,3);
            DDkappa1.row(3).segment(0,3) = - DDkappa1.row(3).segment(0,3);
            DDkappa1.row(3).segment(4,3) = - DDkappa1.row(3).segment(4,3);
            DDkappa1.row(3).segment(8,3) = - DDkappa1.row(3).segment(8,3);

            DDkappa2.col(3).segment(0,3) = - DDkappa2.col(3).segment(0,3);
            DDkappa2.col(3).segment(4,3) = - DDkappa2.col(3).segment(4,3);
            DDkappa2.col(3).segment(8,3) = - DDkappa2.col(3).segment(8,3);
            DDkappa2.row(3).segment(0,3) = - DDkappa2.row(3).segment(0,3);
            DDkappa2.row(3).segment(4,3) = - DDkappa2.row(3).segment(4,3);
            DDkappa2.row(3).segment(8,3) = - DDkappa2.row(3).segment(8,3);
        }

        if ( plate->v_bendingElement[i].sign_2 < 0 )
        {
            DDkappa1.col(7).segment(0,3) = - DDkappa1.col(7).segment(0,3);
            DDkappa1.col(7).segment(4,3) = - DDkappa1.col(7).segment(4,3);
            DDkappa1.col(7).segment(8,3) = - DDkappa1.col(7).segment(8,3);
            DDkappa1.row(7).segment(0,3) = - DDkappa1.row(7).segment(0,3);
            DDkappa1.row(7).segment(4,3) = - DDkappa1.row(7).segment(4,3);
            DDkappa1.row(7).segment(8,3) = - DDkappa1.row(7).segment(8,3);

            DDkappa2.col(7).segment(0,3) = - DDkappa2.col(7).segment(0,3);
            DDkappa2.col(7).segment(4,3) = - DDkappa2.col(7).segment(4,3);
            DDkappa2.col(7).segment(8,3) = - DDkappa2.col(7).segment(8,3);
            DDkappa2.row(7).segment(0,3) = - DDkappa2.row(7).segment(0,3);
            DDkappa2.row(7).segment(4,3) = - DDkappa2.row(7).segment(4,3);
            DDkappa2.row(7).segment(8,3) = - DDkappa2.row(7).segment(8,3);
        }

        len = plate->v_bendingElement[i].voroniLength;
        relevantPart.col(0)=gradKappa1.row(i);
        relevantPart.col(1)=gradKappa2.row(i);

        Jbb = - 1.0/len * relevantPart * EIMat * relevantPart.transpose();

        kappaL = plate->v_bendingElement[i].kappa - plate->v_bendingElement[i].kappaBar;

        temp = - 1.0 / len * kappaL.transpose() * EIMat;
        
        Jbb = Jbb + temp(0) * DDkappa1 + temp(1) * DDkappa2;

        localDOF = plate->v_bendingElement[i].arrayNum;

		for (int j = 0; j < 11; j++)
        {
            for (int k = 0; k < 11; k++)
            {
				int ind1 = localDOF(j);
				int ind2 = localDOF(k);
				stepper->addJacobian(ind1, ind2, - Jbb(k,j));
            }
        }
    }
    
}

// Utility
void elasticBendingForce::crossMat(const Vector3d &a,Matrix3d &b)
{
	b<<0,-a(2),a(1),
	a(2),0,-a(0),
	-a(1),a(0),0;
}

void elasticBendingForce::setFirstJacobian()
{
    for (int k = 0; k < plate->bendingNum; k++)
    {
        localDOF = plate->v_bendingElement[k].arrayNum;

        for (int i = 0; i < 11; i++)
        {
            for (int j = 0; j < 11; j++)
            {
                stepper->addJacobian(localDOF(i), localDOF(j), 1);
            }
        }
    }
}