#ifndef _INDEX_DEFINES_H_
#define _INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>
#include "mechanism.h"

#define URHO 0
#define UMX 1
#define UMY 2
#define UMZ 3
#define UEDEN 4
#define UEINT 5
#define UTEMP 6
#define UFA 7

#define QRHO 0
#define QU 1
#define QV 2
#define QW 3
#define QGAME 4
#define QPRES 5
#define QREINT 6
#define QTEMP 7
#define QFA 8

#define QGAMC 0
#define QC 1
#define QCSML 2
#define QDPDR 3
#define QDPDE 4
#define QRSPEC 5

#define GDRHO 0
#define GDU 1
#define GDV 2
#define GDW 3
#define GDPRES 4
#define GDGAME 5

#define dComp_rhoD 0
#define dComp_mu dComp_rhoD + NUM_SPECIES
#define dComp_xi dComp_mu + 1
#define dComp_lambda dComp_xi + 1

#define comp_Cs2 0
#define comp_CI comp_Cs2 + 1
#define comp_PrT comp_CI + 1
#define comp_Cs2ovPrT comp_PrT
#define nCompC comp_PrT + 1

#define SlipWall 4
#define NoSlipWall 5

#define NTHERM 7
#define QTHERM 8

#ifndef PELEC_USE_AUX
#define NUM_ADV 0
#define NUM_AUX 0
#endif

#define UFS UFA + NUM_ADV
#define UFX UFS + NUM_SPECIES
#define NVAR NTHERM + NUM_ADV + NUM_SPECIES + NUM_AUX

#define QFS QFA + NUM_ADV
#define QFX QFS + NUM_SPECIES
#define QVAR QTHERM + NUM_ADV + NUM_SPECIES + NUM_AUX
#define NQ QVAR
#define NQAUX 6
#define NGDNV 6

#define NPASSIVE NUM_SPECIES + NUM_ADV + NUM_AUX + 3 - AMREX_SPACEDIM

// This value was hardcoded to 4 in get_method_params in the old Fortran code
#define NUM_GROW 4

namespace indxmap {

extern AMREX_GPU_DEVICE_MANAGED int qpassMap[NPASSIVE];
extern AMREX_GPU_DEVICE_MANAGED int upassMap[NPASSIVE];

void init();

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
qpass_map(const int i)
{
  return qpassMap[i];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
upass_map(const int i)
{
  return upassMap[i];
}
} // namespace indxmap

#endif
