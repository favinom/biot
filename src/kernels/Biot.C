//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Biot.h"

registerMooseObject("biotApp", Biot);

defineLegacyParams(Biot);

InputParameters
Biot::validParams()
{
  InputParameters params = Kernel::validParams();
  //params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
  //                           "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  params.addRequiredCoupledVar("disp_x", "real part of the first  coupled component");
  params.addCoupledVar        ("disp_y", "real part of the second coupled component");
  params.addCoupledVar        ("disp_z", "real part of the third  coupled component");
  params.addRequiredCoupledVar("pres"  , "real part of pressure");

  return params;
}

Biot::Biot(const InputParameters & parameters) :
Kernel(parameters),
_dim(_mesh.dimension()),
_disp_x_var(            coupled("disp_x")         ),
_disp_y_var( _dim > 1 ? coupled("disp_y") : 100000),
_disp_z_var( _dim ==3 ? coupled("disp_z") : 100000),
_p_var(coupled("pres")),
_mu(getMaterialProperty<Real>("_mu")),
_k(getMaterialProperty<Real>("_k")),
_kappa(getMaterialProperty<Real>("_kappa")),
_trE(getMaterialProperty<Real>("_divU")),
_trEold(getMaterialProperty<Real>("_divUold")),
_sigma(getMaterialProperty<RealTensorValue>("_sigma"))
{
	_local_i_var=getLocalIndex( _var.number() );

	_id=RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	for (int j=0; j<_dim; ++j)
	{
		_id(j,j)=1.0;
	}
}

Real
Biot::computeQpResidual()
{
	if ( 0<=_local_i_var && _local_i_var<=2)
	{
		initTensor(_V, _grad_test[_i][_qp], _local_i_var);
		_returnValue=_sigma[_qp].contract(_V);
	}
	else
	{
		_returnValue= -(_trE[_qp]-_trEold[_qp]) * _test[_i][_qp];
		if (_is_transient)
		{
			_returnValue-=_dt*_kappa[_qp]*_grad_u[_qp] * _grad_test[_i][_qp];
		}
		else
		{
			_returnValue-= _kappa[_qp]*_grad_u[_qp] * _grad_test[_i][_qp];	
		}
	}
	return _returnValue;
}

Real
Biot::computeQpJacobian()
{
	if ( 0<=_local_i_var && _local_i_var<=2)
	{
		initTensor(_V, _grad_test[_i][_qp], _local_i_var);
		initTensor(_H, _grad_phi [_j][_qp], _local_i_var);

		Real trV=_V.tr();
		Real trH=_H.tr();
		
		RealTensorValue DH=_H-(trH/_dim)*_id;
		RealTensorValue DV=_V-(trV/_dim)*_id;
        
        _returnValue=2.0*_mu[_qp]*DH.contract(DV)+_k[_qp]*trV*trH;
	}
	else
	{
		if (_is_transient)
		{
			_returnValue= -_dt*_kappa[_qp]*_grad_phi[_j][_qp] * _grad_test[_i][_qp];
		}
		else
		{
			_returnValue= -    _kappa[_qp]*_grad_phi[_j][_qp] * _grad_test[_i][_qp];
		}
	}
	return _returnValue;
}


Real Biot::computeQpOffDiagJacobian(unsigned int jvar)
{
	_local_j_var=getLocalIndex(jvar);

	if ( 0<=_local_i_var && _local_i_var<=2)
	{
		initTensor(_V, _grad_test[_i][_qp], _local_i_var);
		Real trV=_V.tr();
		RealTensorValue DV=_V-(trV/_dim)*_id;

		if ( 0<=_local_j_var && _local_j_var<=2)
		{
			initTensor(_H, _grad_phi [_j][_qp], _local_j_var);
			Real trH=_H.tr();
			RealTensorValue DH=_H-(trH/_dim)*_id;			
			
			_returnValue=2.0*_mu[_qp]*DH.contract(DV)+_k[_qp]*trV*trH;
		}
		else
		{
			_returnValue=-_phi[_j][_qp]*trV;
		}
        
	}
	else
	{
		if ( 0<=_local_j_var && _local_j_var<=2)
		{
			initTensor(_H, _grad_phi [_j][_qp], _local_j_var);
			Real trH=_H.tr();
			_returnValue=-trH*_test[_i][_qp];
		}
	}
	return _returnValue;
}

unsigned int Biot::getLocalIndex(unsigned int global_var) const
{

	int _local_var=4;

	if (global_var == _disp_x_var)
		_local_var=0;
	if (global_var == _disp_y_var)
		_local_var=1;
	if (global_var == _disp_z_var)
		_local_var=2;

	if (global_var == _p_var)
		_local_var=3;

	if (_local_var==4)
	{
		std::cout<<"_local_var=4, exiting...\n";
		exit(1);
	}
	return _local_var;
}

void Biot::initTensor(RealTensorValue & V, RealVectorValue const & grad, unsigned int row)
{
		V=RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
		for (int j=0; j<_dim; ++j)
		{
			V(row,j)=grad(j);
		}
		V=0.5*( V+V.transpose() );	
}