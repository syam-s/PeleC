#include "IndexDefines.H"

namespace indxmap {
//extern amrex::GpuArray<int, NPASSIVE> upassMap;
//extern amrex::GpuArray<int, NPASSIVE> qpassMap;

void
init()
{
/*  AMREX_D_PICK(upassMap[0] = UMY; qpassMap[0] = QV; upassMap[1] = UMZ;
               qpassMap[1] = QW; int curMapIndx = 2;, upassMap[0] = UMZ;
               qpassMap[0] = QW; int curMapIndx = 1;, int curMapIndx = 0;);
  for (int i = 0; i < NUM_ADV; ++i) {
    upassMap[curMapIndx] = i + UFA;
    qpassMap[curMapIndx] = i + QFA;
    curMapIndx++;
  }
  for (int i = 0; i < NUM_SPECIES; ++i) {
    upassMap[curMapIndx] = i + UFS;
    qpassMap[curMapIndx] = i + QFS;
    curMapIndx++;
  }*/
}
} // namespace indxmap
