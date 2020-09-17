#ifndef _TIMESTEP_H_
#define _TIMESTEP_H_

#include <AMReX_FArrayBox.H>
#ifdef PELEC_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#endif

#include "EOS.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "Transport.H"

// EstDt routines

namespace TimeStep {
constexpr amrex::Real max_dt = 1.e37;
} // namespace TimeStep

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pc_trans4dt(
  const int which_trans,
  amrex::Real T,
  amrex::Real rho,
  amrex::Real massfrac[],
  amrex::Real& D)
{
  bool get_xi = false, get_mu = false, get_lam = false, get_Ddiag = false;
  amrex::Real dum1 = 0., dum2 = 0.;

  if (which_trans == 0) {
    get_mu = true;
    transport::transport(
      get_xi, get_mu, get_lam, get_Ddiag, T, rho, massfrac, nullptr, D, dum1,
      dum2);
  } else if (which_trans == 1) {
    get_lam = true;
    transport::transport(
      get_xi, get_mu, get_lam, get_Ddiag, T, rho, massfrac, nullptr, dum1, dum2,
      D);
  }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
pc_estdt_hydro(
  amrex::Box const& bx,
  const amrex::Array4<const amrex::Real>& u,
#ifdef PELEC_USE_EB
  const amrex::Array4<const amrex::EBCellFlag>& flags,
#endif
  AMREX_D_DECL(
    const amrex::Real& dx,
    const amrex::Real& dy,
    const amrex::Real& dz)) noexcept
{
  amrex::Real dt = TimeStep::max_dt;

  amrex::Loop(bx, [=, &dt](int i, int j, int k) {
#ifdef PELEC_USE_EB
    if (not flags(i, j, k).isCovered()) {
#endif
      const amrex::Real rho = u(i, j, k, URHO);
      const amrex::Real rhoInv = 1.0 / rho;
      amrex::Real T = u(i, j, k, UTEMP);
      amrex::Real massfrac[NUM_SPECIES];
      amrex::Real c;
      for (int n = 0; n < NUM_SPECIES; ++n)
        massfrac[n] = u(i, j, k, UFS + n) * rhoInv;
      EOS::RTY2Cs(rho, T, massfrac, c);
      AMREX_D_TERM(const amrex::Real ux = u(i, j, k, UMX) * rhoInv;
                   const amrex::Real dt1 = dx / (c + amrex::Math::abs(ux));
                   dt = amrex::min(dt, dt1);
                   , const amrex::Real uy = u(i, j, k, UMY) * rhoInv;
                   const amrex::Real dt2 = dy / (c + amrex::Math::abs(uy));
                   dt = amrex::min(dt, dt2);
                   , const amrex::Real uz = u(i, j, k, UMZ) * rhoInv;
                   const amrex::Real dt3 = dz / (c + amrex::Math::abs(uz));
                   dt = amrex::min(dt, dt3););
#ifdef PELEC_USE_EB
    }
#endif
  });
  return dt;
}

// Diffusion Velocity
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
pc_estdt_veldif(
  amrex::Box const bx,
  const amrex::Array4<const amrex::Real>& u,
#ifdef PELEC_USE_EB
  const amrex::Array4<const amrex::EBCellFlag>& flags,
#endif
  AMREX_D_DECL(
    const amrex::Real& dx,
    const amrex::Real& dy,
    const amrex::Real& dz)) noexcept
{
  amrex::Real dt = TimeStep::max_dt;
  int which_trans = 0;

  amrex::Loop(bx, [=, &dt](int i, int j, int k) {
#ifdef PELEC_USE_EB
    if (not flags(i, j, k).isCovered()) {
#endif
      const amrex::Real rho = u(i, j, k, URHO);
      const amrex::Real rhoInv = 1.0 / rho;
      amrex::Real massfrac[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n) {
        massfrac[n] = u(i, j, k, n + UFS) * rhoInv;
      }
      amrex::Real T = u(i, j, k, UTEMP);
      amrex::Real D = 0.0;
      pc_trans4dt(which_trans, T, rho, massfrac, D);
      D *= rhoInv;
      if (D == 0.0)
        D = SMALL;
      AMREX_D_TERM(
        const amrex::Real dt1 = 0.5 * dx * dx / (AMREX_SPACEDIM * D);
        dt = amrex::min(dt, dt1);
        , const amrex::Real dt2 = 0.5 * dy * dy / (AMREX_SPACEDIM * D);
        dt = amrex::min(dt, dt2);
        , const amrex::Real dt3 = 0.5 * dz * dz / (AMREX_SPACEDIM * D);
        dt = amrex::min(dt, dt3););
#ifdef PELEC_USE_EB
    }
#endif
  });
  return dt;
}

// Diffusion Temperature
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
pc_estdt_tempdif(
  amrex::Box const bx,
  const amrex::Array4<const amrex::Real>& u,
#ifdef PELEC_USE_EB
  const amrex::Array4<const amrex::EBCellFlag>& flags,
#endif
  AMREX_D_DECL(
    const amrex::Real& dx,
    const amrex::Real& dy,
    const amrex::Real& dz)) noexcept
{
  amrex::Real dt = TimeStep::max_dt;
  int which_trans = 1;

  amrex::Loop(bx, [=, &dt](int i, int j, int k) {
#ifdef PELEC_USE_EB
    if (not flags(i, j, k).isCovered()) {
#endif
      const amrex::Real rho = u(i, j, k, URHO);
      const amrex::Real rhoInv = 1.0 / rho;
      amrex::Real massfrac[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n)
        massfrac[n] = u(i, j, k, n + UFS) * rhoInv;
      amrex::Real T = u(i, j, k, UTEMP);
      amrex::Real D = 0.0;
      pc_trans4dt(which_trans, T, rho, massfrac, D);
      amrex::Real cv;
      EOS::TY2Cv(T, massfrac, cv);
      D *= rhoInv / cv;
      if (D == 0.0)
        D = SMALL;
      AMREX_D_TERM(
        const amrex::Real dt1 = 0.5 * dx * dx / (AMREX_SPACEDIM * D);
        dt = amrex::min(dt, dt1);
        , const amrex::Real dt2 = 0.5 * dy * dy / (AMREX_SPACEDIM * D);
        dt = amrex::min(dt, dt2);
        , const amrex::Real dt3 = 0.5 * dz * dz / (AMREX_SPACEDIM * D);
        dt = amrex::min(dt, dt3););
#ifdef PELEC_USE_EB
    }
#endif
  });
  return dt;
}

// Diffusion Enthalpy
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
pc_estdt_enthdif(
  amrex::Box const bx,
  const amrex::Array4<const amrex::Real>& u,
#ifdef PELEC_USE_EB
  const amrex::Array4<const amrex::EBCellFlag>& flags,
#endif
  AMREX_D_DECL(
    const amrex::Real& dx,
    const amrex::Real& dy,
    const amrex::Real& dz)) noexcept
{
  amrex::Real dt = TimeStep::max_dt;
  int which_trans = 1;

  amrex::Loop(bx, [=, &dt](int i, int j, int k) {
#ifdef PELEC_USE_EB
    if (not flags(i, j, k).isCovered()) {
#endif
      const amrex::Real rho = u(i, j, k, URHO);
      const amrex::Real rhoInv = 1.0 / rho;
      amrex::Real massfrac[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n)
        massfrac[n] = u(i, j, k, n + UFS) * rhoInv;
      amrex::Real T = u(i, j, k, UTEMP);
      amrex::Real cp;
      EOS::TY2Cp(T, massfrac, cp);
      amrex::Real D;
      pc_trans4dt(which_trans, T, rho, massfrac, D);
      D *= rhoInv / cp;
      AMREX_D_TERM(
        const amrex::Real dt1 = 0.5 * dx * dx / (AMREX_SPACEDIM * D);
        dt = amrex::min(dt, dt1);
        , const amrex::Real dt2 = 0.5 * dy * dy / (AMREX_SPACEDIM * D);
        dt = amrex::min(dt, dt2);
        , const amrex::Real dt3 = 0.5 * dz * dz / (AMREX_SPACEDIM * D);
        dt = amrex::min(dt, dt3););

#ifdef PELEC_USE_EB
    }
#endif
  });
  return dt;
}

#endif
