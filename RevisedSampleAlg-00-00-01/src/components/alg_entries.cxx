#include "GaudiKernel/DeclareFactoryEntries.h"
#include "SampleAlg/JpsiToPhiEtaAlg.h"

#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/SecondVertexFit.h"

VertexFit * vtxfit = VertexFit::instance();
KinematicFit * kmfit = KinematicFit::instance();
SecondVertexFit *ksfit = SecondVertexFit::instance();

DECLARE_ALGORITHM_FACTORY( JpsiToKsKlEtaAlg )

DECLARE_FACTORY_ENTRIES( RevisedSampleAlg ) {
  DECLARE_ALGORITHM(JpsiToKsKlEtaAlg)
}


