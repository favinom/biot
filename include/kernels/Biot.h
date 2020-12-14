//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

class Biot;

template <>
InputParameters validParams<Biot>();

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class Biot : public Kernel
{
public:
  static InputParameters validParams();

  Biot(const InputParameters & parameters);

  unsigned int getLocalIndex(unsigned int) const;
  void initTensor(RealTensorValue & V, RealVectorValue const & grad, unsigned int row);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  virtual Real computeQpOffDiagJacobian(unsigned int /*jvar*/) override;

  unsigned int _dim;
  unsigned int _disp_x_var;
  unsigned int _disp_y_var;
  unsigned int _disp_z_var;
  unsigned int _p_var;

  MaterialProperty<Real>            const & _mu;
  MaterialProperty<Real>            const & _k;
  MaterialProperty<Real>            const & _kappa;
  MaterialProperty<Real>            const & _trE;
  MaterialProperty<Real>            const & _trEold;
  MaterialProperty<RealTensorValue> const & _sigma;

  RealTensorValue _id;

  unsigned int _local_i_var;
  unsigned int _local_j_var;
  Real _returnValue;
  RealTensorValue _V;
  RealTensorValue _H;

};
