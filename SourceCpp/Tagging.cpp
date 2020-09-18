#include <AMReX_ParmParse.H>

#include "PeleC.H"
#include "Tagging.H"

void
PeleC::read_tagging_params(
  amrex::GpuArray<amrex::Real, max_prob_param>& tagging_parm_real,
  amrex::GpuArray<int, max_prob_param>& tagging_parm_int)
{
  tagging_parm_real[denerr] = 1.0e10;
  tagging_parm_real[presserr] = 1.0e10;
  tagging_parm_real[velerr] = 1.0e10;
  tagging_parm_real[vorterr] = 1.0e10;
  tagging_parm_real[temperr] = 1.0e10;
  tagging_parm_real[ftracerr] = 1.0e10;
  tagging_parm_real[vfracerr] = 1.0e10;
  tagging_parm_real[dengrad] = 1.0e10;
  tagging_parm_real[pressgrad] = 1.0e10;
  tagging_parm_real[velgrad] = 1.0e10;
  tagging_parm_real[tempgrad] = 1.0e10;
  tagging_parm_real[ftracgrad] = 1.0e10;

  tagging_parm_int[max_denerr_lev] = 10;
  tagging_parm_int[max_presserr_lev] = 10;
  tagging_parm_int[max_velerr_lev] = 10;
  tagging_parm_int[max_vorterr_lev] = 10;
  tagging_parm_int[max_temperr_lev] = 10;
  tagging_parm_int[max_ftracerr_lev] = 10;
  tagging_parm_int[max_vfracerr_lev] = 10;
  tagging_parm_int[max_dengrad_lev] = 10;
  tagging_parm_int[max_pressgrad_lev] = 10;
  tagging_parm_int[max_velgrad_lev] = 10;
  tagging_parm_int[max_tempgrad_lev] = 10;
  tagging_parm_int[max_ftracgrad_lev] = 10;

  amrex::ParmParse pp("tagging");

  pp.query("denerr", tagging_parm_real[denerr]);
  pp.query("presserr", tagging_parm_real[presserr]);
  pp.query("velerr", tagging_parm_real[velerr]);
  pp.query("vorterr", tagging_parm_real[vorterr]);
  pp.query("temperr", tagging_parm_real[temperr]);
  pp.query("ftracerr", tagging_parm_real[ftracerr]);
  pp.query("vfracerr", tagging_parm_real[vfracerr]);

  pp.query("max_denerr_lev", tagging_parm_int[max_denerr_lev]);
  pp.query("max_presserr_lev", tagging_parm_int[max_presserr_lev]);
  pp.query("max_velerr_lev", tagging_parm_int[max_velerr_lev]);
  pp.query("max_vorterr_lev", tagging_parm_int[max_vorterr_lev]);
  pp.query("max_temperr_lev", tagging_parm_int[max_temperr_lev]);
  pp.query("max_ftracerr_lev", tagging_parm_int[max_ftracerr_lev]);
  pp.query("max_vfracerr_lev", tagging_parm_int[max_vfracerr_lev]);

  pp.query("dengrad", tagging_parm_real[dengrad]);
  pp.query("pressgrad", tagging_parm_real[pressgrad]);
  pp.query("velgrad", tagging_parm_real[velgrad]);
  pp.query("tempgrad", tagging_parm_real[tempgrad]);
  pp.query("ftracgrad", tagging_parm_real[ftracgrad]);

  pp.query("max_dengrad_lev", tagging_parm_int[max_dengrad_lev]);
  pp.query("max_pressgrad_lev", tagging_parm_int[max_pressgrad_lev]);
  pp.query("max_velgrad_lev", tagging_parm_int[max_velgrad_lev]);
  pp.query("max_tempgrad_lev", tagging_parm_int[max_tempgrad_lev]);
  pp.query("max_ftracgrad_lev", tagging_parm_int[max_ftracgrad_lev]);
}
