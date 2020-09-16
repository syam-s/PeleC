#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "mechanism.h"

#include "IndexDefines.H"
#include "Constants.H"
#include "EOS.H"
#include "Tagging.H"
#include "Transport.H"
#include "ProblemDerive.H"
#include "PeleC.H"

enum ProbCompReal {
  reynolds = 0,
  mach,
  prandtl,
  omega_x,
  omega_y,
  omega_z,
  L_x,
  L_y,
  L_z,
  L,
  p0,
  T0,
  rho0,
  v0
};

enum ProbCompBool {
  convecting = 0
};

static void prob_param_fill(amrex::GeometryData const& geomdata, amrex::GpuArray<amrex::Real,max_prob_param>& prob_parm_real, amrex::GpuArray<bool,max_prob_param>& prob_parm_bool)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();

  prob_parm_real[reynolds] = 1600.0;
  prob_parm_real[mach] = 0.1;
  prob_parm_real[prandtl] = 0.71;
  prob_parm_real[omega_x] = 1.0;
  prob_parm_real[omega_y] = 1.0;
  prob_parm_real[omega_z] = 1.0;
  prob_parm_real[L_x] = 0.0;
  prob_parm_real[L_y] = 0.0;
  prob_parm_real[L_z] = 0.0;
  prob_parm_real[L] = 0.0;
  prob_parm_real[p0] = 1.013e6;
  prob_parm_real[T0] = 300.0;
  prob_parm_real[rho0] = 0.0;
  prob_parm_real[v0] = 0.0;
  prob_parm_bool[convecting] = false;

  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("reynolds", prob_parm_real[reynolds]);
  pp.query("mach", prob_parm_real[mach]);
  pp.query("prandtl", prob_parm_real[prandtl]);
  pp.query("omega_x", prob_parm_real[omega_x]);
  pp.query("omega_y", prob_parm_real[omega_y]);
  pp.query("omega_z", prob_parm_real[omega_z]);
  pp.query("convecting", prob_parm_bool[convecting]);

  // Define the length scale
  prob_parm_real[L] = 1.0 / PI;
  prob_parm_real[L_x] = prob_hi[0] - prob_lo[0];
  prob_parm_real[L_y] = prob_hi[1] - prob_lo[1];
  prob_parm_real[L_z] = prob_hi[2] - prob_lo[2];

  // Initial density, velocity, and material properties
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  EOS::PYT2RE(prob_parm_real[p0], massfrac, prob_parm_real[T0], prob_parm_real[rho0], eint);
  EOS::RTY2Cs(prob_parm_real[rho0], prob_parm_real[T0], massfrac, cs);
  EOS::TY2Cp(prob_parm_real[T0], massfrac, cp);

  prob_parm_real[v0] = prob_parm_real[mach] * cs;
  transport_params::const_bulk_viscosity = 0.0;
  transport_params::const_diffusivity = 0.0;
  transport_params::const_viscosity =
    prob_parm_real[rho0] * prob_parm_real[v0] * prob_parm_real[L] / prob_parm_real[reynolds];
  transport_params::const_conductivity =
    transport_params::const_viscosity * cp / prob_parm_real[prandtl];

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "L, rho0, v0, p0, T0, gamma, mu, k, c_s0, Reynolds, "
                       "Mach, Prandtl, omega_x, omega_y, omega_z"
                    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << prob_parm_real[L] << "," << prob_parm_real[rho0] << "," << prob_parm_real[v0] << ","
    << prob_parm_real[p0] << "," << prob_parm_real[T0] << "," << EOS::gamma << ","
    << transport_params::const_viscosity << ","
    << transport_params::const_conductivity << "," << cs << ","
    << prob_parm_real[reynolds] << "," << prob_parm_real[mach] << "," << prob_parm_real[prandtl]
    << "," << prob_parm_real[omega_x] << "," << prob_parm_real[omega_y] << ","
    << prob_parm_real[omega_z] << std::endl;
  ofs.close();
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  const amrex::GpuArray<amrex::Real,max_prob_param>& prob_parm_real,
  const amrex::GpuArray<bool,max_prob_param>& prob_parm_bool)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  amrex::Real z = 0;
#if AMREX_SPACEDIM == 3
  z = prob_lo[2] + (k + 0.5) * dx[2];
#endif

  // TG functions
  amrex::Real u[3] = {0.0};
  u[0] = prob_parm_real[v0] * sin(prob_parm_real[omega_x] * x / prob_parm_real[L]) *
         cos(prob_parm_real[omega_y] * y / prob_parm_real[L]) *
         cos(prob_parm_real[omega_z] * z / prob_parm_real[L]);
  u[1] = -prob_parm_real[v0] * cos(prob_parm_real[omega_x] * x / prob_parm_real[L]) *
         sin(prob_parm_real[omega_y] * y / prob_parm_real[L]) *
         cos(prob_parm_real[omega_z] * z / prob_parm_real[L]);
  if (prob_parm_bool[convecting]) {
    u[0] += prob_parm_real[v0];
    u[1] += prob_parm_real[v0];
  }
  const amrex::Real p =
    prob_parm_real[p0] + prob_parm_real[rho0] * prob_parm_real[v0] * prob_parm_real[v0] / 16.0 *
                     (cos(2.0 * prob_parm_real[omega_x] * x / prob_parm_real[L]) +
                      cos(2.0 * prob_parm_real[omega_y] * y / prob_parm_real[L])) *
                     (cos(2.0 * prob_parm_real[omega_z] * z / prob_parm_real[L]) + 2.0);
  amrex::Real rho, eint;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  EOS::PYT2RE(p, massfrac, prob_parm_real[T0], rho, eint);

  // Set the state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u[0];
  state(i, j, k, UMY) = rho * u[1];
  state(i, j, k, UMZ) = rho * u[2];
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) =
    rho * (eint + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
  state(i, j, k, UTEMP) = prob_parm_real[T0];
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = rho * massfrac[n];
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata)
{
}

struct MyProbTagStruct
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void set_problem_tags(
    const int i,
    const int j,
    const int k,
    amrex::Array4<char> const& tag,
    amrex::Array4<amrex::Real const> const& field,
    char tagval,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo,
    const amrex::Real time,
    const int level) noexcept
  {
    // could do problem specific tagging here
  }
};

using ProblemTags = MyProbTagStruct;

struct MyProbDeriveStruct
{
  static void
  add(amrex::DeriveList& derive_lst, amrex::DescriptorList& desc_lst)
  {
    // Add derives as follows and define the derive function below:
    // derive_lst.add(
    //  "varname", amrex::IndexType::TheCellType(), 1, pc_varname,
    //  the_same_box);
    // derive_lst.addComponent("varname", desc_lst, State_Type, 0, NVAR);
  }

  static void pc_varname(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int /*ncomp*/,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real /*time*/,
    const int* /*bcrec*/,
    int level)
  {
    auto const dat = datfab.array();
    auto arr = derfab.array();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // do something with arr
    });
  }
};

void pc_prob_close();

using ProblemDerives = MyProbDeriveStruct;

#endif
