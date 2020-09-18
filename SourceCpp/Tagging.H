#ifndef _TAGGING_H_
#define _TAGGING_H_

#include <cmath>

#include <AMReX_FArrayBox.H>
#include <AMReX_TagBox.H>

enum TaggingCompReal {
  denerr = 0,
  presserr,
  velerr,
  vorterr,
  temperr,
  ftracerr,
  vfracerr,
  dengrad,
  pressgrad,
  velgrad,
  tempgrad,
  ftracgrad
};

enum TaggingCompInt {
  max_denerr_lev = 0,
  max_presserr_lev,
  max_velerr_lev,
  max_temperr_lev,
  max_ftracerr_lev,
  max_vfracerr_lev,
  max_dengrad_lev,
  max_pressgrad_lev,
  max_velgrad_lev,
  max_vorterr_lev,
  max_tempgrad_lev,
  max_ftracgrad_lev
};

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
tag_error(
  const int i,
  const int j,
  const int k,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fielderr,
  char tagval) noexcept
{
  // Tag on regions of high field values
  if (field(i, j, k) >= fielderr) {
    tag(i, j, k) = tagval;
  }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
tag_graderror(
  const int i,
  const int j,
  const int k,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fieldgrad,
  char tagval) noexcept
{
  // Tag on regions of high field gradient
  AMREX_D_TERM(
    amrex::Real ax = amrex::Math::abs(field(i + 1, j, k) - field(i, j, k));
    ax = amrex::max(ax, amrex::Math::abs(field(i, j, k) - field(i - 1, j, k)));
    , amrex::Real ay = amrex::Math::abs(field(i, j + 1, k) - field(i, j, k));
    ay = amrex::max(ay, amrex::Math::abs(field(i, j, k) - field(i, j - 1, k)));
    amrex::Real az = amrex::Math::abs(field(i, j, k + 1) - field(i, j, k));
    ,
    az = amrex::max(az, amrex::Math::abs(field(i, j, k) - field(i, j, k - 1)));)
#if AMREX_SPACEDIM > 1
  if (amrex::max(AMREX_D_DECL(ax, ay, az)) >= fieldgrad)
#elif AMREX_SPACEDIM == 1
  if (ax >= fieldgrad)
#endif
  {
    tag(i, j, k) = tagval;
  }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
tag_abserror(
  const int i,
  const int j,
  const int k,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fielderr,
  char tagval) noexcept
{
  // Tag on regions of high field values
  if (amrex::Math::abs(field(i, j, k)) >= fielderr) {
    tag(i, j, k) = tagval;
  }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
tag_error_bounds(
  const int i,
  const int j,
  const int k,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real lbnd,
  const amrex::Real ubnd,
  char tagval) noexcept
{
  // Tag on regions inside bounds
  if ((lbnd < field(i, j, k)) and (field(i, j, k) < ubnd)) {
    tag(i, j, k) = tagval;
  }
}

struct EmptyProbTagStruct
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
  }
};

template <typename ProbTagStruct, typename... Args>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
set_problem_tags(Args&&... args)
{
  ProbTagStruct::set_problem_tags(std::forward<Args>(args)...);
}

#endif
