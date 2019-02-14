#include "stdafx.h"

void GetEigenVector(const PNEANet& Graph, TIntFltH& NIdEigenH, const double& Eps=1e-4, const int& MaxIter=100) {
  const int NNodes = Graph->GetNodes();
  NIdEigenH.Gen(NNodes);
  // initialize vector values
  for (TNEANet::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    NIdEigenH.AddDat(NI.GetId(), 1.0/NNodes);
    IAssert(NI.GetId() == NIdEigenH.GetKey(NIdEigenH.Len()-1));
  }
  TFltV TmpV(NNodes);
  for (int iter = 0; iter < MaxIter; iter++) {
    int j = 0;
    // add neighbor values
    for (TNEANet::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++, j++) {
      TmpV[j] = 0; const int NodeId = NI.GetId();
      for (int e = 0; e < NI.GetOutDeg(); e++) {
        const int NeigNId = NI.GetOutNId(e);
        const int EId = Graph->GetEId(NodeId, NeigNId);
        const double weight = Graph->GetFltAttrDatE(EId, TSnap::CapAttrName);
        TmpV[j] += weight + NIdEigenH.GetDat(NeigNId); }
      for (int e = 0; e < NI.GetInDeg(); e++) {
        const int NeigNId = NI.GetInNId(e);
        const int EId = Graph->GetEId(NeigNId, NodeId);
        const double weight = Graph->GetFltAttrDatE(EId, TSnap::CapAttrName);
        TmpV[j] += weight + NIdEigenH.GetDat(NeigNId); }
    }

    // normalize
    double sum = 0;
    for (int i = 0; i < TmpV.Len(); i++) {
      sum += (TmpV[i]*TmpV[i]);
    }
    sum = sqrt(sum);
    for (int i = 0; i < TmpV.Len(); i++) {
      TmpV[i] /= sum;
    }

    // compute difference
    double diff = 0.0;
    j = 0;
    for (TNEANet::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++, j++) {
      diff += fabs(NIdEigenH.GetDat(NI.GetId())-TmpV[j]);
    }

    // set new values
    j = 0;
    for (TNEANet::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++, j++) {
      NIdEigenH.AddDat(NI.GetId(), TmpV[j]);
    }

    if (diff < Eps) {
      break;
    }
  }
}

//### build a directed network graph with weight value on each edge
void BuildNetwork(const TStr& InFNm, PNEANet& Net) {
  TSsParser Ss(InFNm, ssfWhiteSep, true, true, true);
  Net.Clr(); Net = TNEANet::New();
  int SrcNId, DstNId, EId; double WeightVal;
  while (Ss.Next()) {
    if (! Ss.GetInt(0, SrcNId) || ! Ss.GetInt(1, DstNId)) { continue; }
    Ss.GetFlt(2, WeightVal);
    if (! Net->IsNode(SrcNId)) { 
      Net->AddNode(SrcNId); //### add source node if it does not exist on the net
    }
    if (! Net->IsNode(DstNId)) {
      Net->AddNode(DstNId); //### add destination node if it does not exist on the net
    }
    EId = Net->AddEdge(SrcNId, DstNId); //### connect the two nodes
    Net->AddFltAttrDatE(EId, WeightVal, TSnap::CapAttrName); //### add weight value b/w nodes
  }
  Net->Defrag(); //### defragment the net (compact the memory)
}

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("Node Centrality. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "../as20graph.txt", "Input un/directed graph");
  //const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "node_centrality.tab", "Output file");
  const TStr OutFNm = "centrality_" + InFNm;
  printf("Loading input: <%s> ... Output: <%s>\n", InFNm.CStr(), OutFNm.CStr());
  PNEANet Graph; BuildNetwork(InFNm, Graph);
  printf("nodes:%d  edges:%d\n", Graph->GetNodes(), Graph->GetEdges());
  //TIntFltH BtwH, CloseH;
  TIntFltH EigH;
  printf("Computing ...\n");
  printf(" Eigenvector...\n");           GetEigenVector(Graph, EigH);
  //printf(" Betweenness (SLOW!)...");   TSnap::GetBetweennessCentr(Graph, BtwH, 1);
  //printf(" Closeness (SLOW!)...");
  //for (TNEANet::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
  //  const int NId = NI.GetId();
  //  CloseH.AddDat(NId, TSnap::GetClosenessCentr<PUNGraph>(UGraph, NId, false));
  //}
  printf("\nDONE! saving...");
  FILE *F = fopen(OutFNm.CStr(), "wt");
  fprintf(F,"#Network: %s\n", InFNm.CStr());
  fprintf(F,"#Nodes: %d\tEdges: %d\n", Graph->GetNodes(), Graph->GetEdges());
  fprintf(F,"#NodeId\tEigenVector\n");
  //fprintf(F,"#NodeId\tCloseness\tBetweennes\tEigenVector\n");
  for (TNEANet::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    const int NId = NI.GetId();
    //const double CloCentr = CloseH.GetDat(NId);
    //const double BtwCentr = BtwH.GetDat(NId);
    const double EigCentr = EigH.GetDat(NId);
    //fprintf(F, "%d\t%f\t%f\t%f\n", NId, CloCentr, BtwCentr, EigCentr);
    fprintf(F, "%d\t%0.10f\n", NId, EigCentr);
  }
  fclose(F);
  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
