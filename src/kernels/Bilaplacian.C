//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Bilaplacian.h"

registerMooseObject("biotApp", Bilaplacian);

defineLegacyParams(Bilaplacian);

InputParameters
Bilaplacian::validParams()
{
  InputParameters params = Kernel::validParams();
  //params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
  //                           "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  params.addRequiredCoupledVar("main_var", "real part of the first  coupled component");
  params.addRequiredCoupledVar("aux_var" , "real part of the second coupled component");

  return params;
}

Bilaplacian::Bilaplacian(const InputParameters & parameters) :
Kernel(parameters),
_main_var( coupled("main_var") ),
_aux_var ( coupled("aux_var")  ),
_main( coupledValue("main_var") ),
_q( coupledValue("aux_var") )
{
	_local_var=4;

	if (_var.number() == _main_var)
		_local_var=0;
	if (_var.number() == _aux_var)
		_local_var=1;

	if (_local_var==4)
	{
		std::cout<<"_local_var=4, exiting...\n";
		exit(1);
	}

}

Real
Bilaplacian::computeQpResidual()
{
	if (_local_var==0)
	{
		return -_grad_u[_qp]*_grad_test[_i][_qp]-_q[_qp]*_test[_i][_qp];
	}
	if (_local_var==1)
	{
		return -_grad_u[_qp]*_grad_test[_i][_qp] + 100.0*_main[_qp]*_test[_i][_qp];
	}
	exit(1);
	return 0.0;
}

Real
Bilaplacian::computeQpJacobian()
{
  //return _grad_phi[_j][_qp] * _grad_test[_i][_qp];
	return 0.0;
}
