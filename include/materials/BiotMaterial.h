//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
//#include "LinearInterpolation.h"

class BiotMaterial;

template <>
InputParameters validParams<BiotMaterial>();

class BiotMaterial : public Material
{
public:
  BiotMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  unsigned int const _dim;
  Real const _mu_input;
  Real const _k_input;
  Real const _kappa_input;

  VariableGradient const & _grad_disp_x;
  VariableGradient const & _grad_disp_y;
  VariableGradient const & _grad_disp_z;
  VariableGradient const & _grad_disp_old_x;
  VariableGradient const & _grad_disp_old_y;
  VariableGradient const & _grad_disp_old_z;

  VariableValue    const & _pres;

  MaterialProperty<Real>            & _mu;
  MaterialProperty<Real>            & _k;
  MaterialProperty<Real>            & _kappa;
  MaterialProperty<Real>            & _trE;
  MaterialProperty<Real>            & _trEold;
  MaterialProperty<RealTensorValue> & _sigma;

  RealTensorValue _id;

};

