#include "biotApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
biotApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

biotApp::biotApp(InputParameters parameters) : MooseApp(parameters)
{
  biotApp::registerAll(_factory, _action_factory, _syntax);
}

biotApp::~biotApp() {}

void
biotApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"biotApp"});
  Registry::registerActionsTo(af, {"biotApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
biotApp::registerApps()
{
  registerApp(biotApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
biotApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  biotApp::registerAll(f, af, s);
}
extern "C" void
biotApp__registerApps()
{
  biotApp::registerApps();
}
