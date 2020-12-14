//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Elasticity.h"

registerMooseObject("biotApp", Elasticity);

defineLegacyParams(Elasticity);

InputParameters
Elasticity::validParams()
{
  InputParameters params = Kernel::validParams();
  //params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
  //                           "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  params.addRequiredCoupledVar("disp_x", "real part of the first  coupled component");
  params.addCoupledVar        ("disp_y", "real part of the second coupled component");
  params.addCoupledVar        ("disp_z", "real part of the third  coupled component");

  return params;
}

Elasticity::Elasticity(const InputParameters & parameters) :
Kernel(parameters),
_dim(_mesh.dimension()),
_disp_x_var(            coupled("disp_x")         ),
_disp_y_var( _dim > 1 ? coupled("disp_y") : 100000),
_disp_z_var( _dim ==3 ? coupled("disp_z") : 100000),
_mu(getMaterialProperty<Real>("_mu")),
_k(getMaterialProperty<Real>("_k")),
_sigma(getMaterialProperty<RealTensorValue>("_sigma"))
{
	_local_var=4;

	if (_var.number() == _disp_x_var)
		_local_var=0;
	if (_var.number() == _disp_y_var)
		_local_var=1;
	if (_var.number() == _disp_z_var)
		_local_var=2;

	if (_local_var==4)
	{
		std::cout<<"_local_var=4, exiting...\n";
		exit(1);
	}

}

Real
Elasticity::computeQpResidual()
{
	if ( 0<=_local_var && _local_var<=2)
	{
		RealTensorValue V(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
		for (int j=0; j<_dim; ++j)
		{
			V(_local_var,j)=_grad_test[_i][_qp](j);
		}
		V=0.5*( V+V.transpose() );
		Real retV=_sigma[_qp].contract(V);
		return retV;
	}
	exit(1);
}

Real
Elasticity::computeQpJacobian()
{
  //return _grad_phi[_j][_qp] * _grad_test[_i][_qp];
	return 0.0;
}
