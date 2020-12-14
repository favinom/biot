  //* This file is part of the MOOSE framework
  //* https://www.mooseframework.org
  //*
  //* All rights reserved, see COPYRIGHT for full restrictions
  //* https://github.com/idaholab/moose/blob/master/COPYRIGHT
  //*
  //* Licensed under LGPL 2.1, please see LICENSE for details
  //* https://www.gnu.org/licenses/lgpl-2.1.html

  #include "BiotMaterial.h"

  registerMooseObject("biotApp", BiotMaterial);

  template <>
  InputParameters
  validParams<BiotMaterial>()
  {
    InputParameters params = validParams<Material>();
    params.addRequiredParam<Real>("mu","mu");
    params.addRequiredParam<Real>("k","k");
    params.addRequiredParam<Real>("kappa","kappa");

    params.addRequiredCoupledVar("disp_x", "displacement in x direction");
    params.addCoupledVar        ("disp_y", "real part of the second coupled component");
    params.addCoupledVar        ("disp_z", "real part of the third  coupled component");
    params.addRequiredCoupledVar("pres","pore pressure");


    return params;
  }

  BiotMaterial::BiotMaterial(const InputParameters & parameters) :
  Material(parameters),
  _dim(_mesh.dimension()),
  _mu_input   ( getParam<Real>("mu")    ),
  _k_input    ( getParam<Real>("k")     ),
  _kappa_input( getParam<Real>("kappa") ),
  _grad_disp_x( coupledGradient("disp_x") ),
  _grad_disp_y( _dim > 1 ? coupledGradient("disp_y") : _grad_zero ),
  _grad_disp_z( _dim ==3 ? coupledGradient("disp_z") : _grad_zero ),
  _grad_disp_old_x( _is_transient             ? coupledGradientOld("disp_x") : _grad_zero ),
  _grad_disp_old_y( _is_transient && _dim > 1 ? coupledGradientOld("disp_y") : _grad_zero ),
  _grad_disp_old_z( _is_transient && _dim ==3 ? coupledGradientOld("disp_z") : _grad_zero ),
  _pres(coupledValue("pres")),
  _mu(declareProperty<Real>("_mu")),
  _k(declareProperty<Real>("_k")),
  _kappa(declareProperty<Real>("_kappa")),
  _trE(declareProperty<Real>("_divU")),
  _trEold(declareProperty<Real>("_divUold")),
  _sigma(declareProperty<RealTensorValue>("_sigma"))
  {
    _id=RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    for (int i=0; i<_dim; ++i)
      _id(i,i)=1.0;
  }

  void BiotMaterial::computeQpProperties()
  {
    _mu[_qp]   =   _mu_input;
    _k[_qp]    =    _k_input;
    _kappa[_qp]=_kappa_input;

    RealTensorValue U(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    RealTensorValue Uold(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

     for (int j=0; j<_dim; ++j)
     {
       U(0,j)=_grad_disp_x[_qp](j);
       U(1,j)=_grad_disp_y[_qp](j);
       U(2,j)=_grad_disp_z[_qp](j);
       Uold(0,j)=_grad_disp_old_x[_qp](j);
       Uold(1,j)=_grad_disp_old_y[_qp](j);
       Uold(2,j)=_grad_disp_old_z[_qp](j);
     }

     RealTensorValue E=0.5*( U+U.transpose() );
     Real trE=E.tr();
     _trE[_qp]=trE;
     _trEold[_qp]=Uold.tr();
     RealTensorValue D=E-(trE/_dim)*_id;

     _sigma[_qp]=2.0*_mu[_qp]*  D + (_k[_qp] * trE -_pres[_qp] )*  _id ;

  }
