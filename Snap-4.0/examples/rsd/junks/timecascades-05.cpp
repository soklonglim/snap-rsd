#include "stdafx.h"
#include <omp.h>
#include <time.h>
#include <math.h>

#define CONSIDER_ORDER

#define DB 51
#define DBR 1
// #define Q_SZ 40
#define Q_Prob 1
#define SIMULATIONSZ 100
// #define InfSrcNode 23 // source of infected node 
#define K_PICK 5
// idev --cpus-per-task=4 --partition=pddms --time=2:00:00
/* Cit-HepTh, CollegeMsg, email-Eu-core, Wiki-Vote, Email-Enron, soc-Epinions */
const TStr INPUTFILE = "dummy.txt";
const TStr OUTPUTFILE = "output-" + INPUTFILE;
const TStr INFECTEDFILE = "si-" + INPUTFILE;

int Q_SZ = 0;

TRnd Rnd(0);
void BuildNetwork(const TStr& InFNm, PNEANet& Net, PNEANet& ReverseNet); //### build directed network graph
void BuildNetworkRandomWeight(const TStr& InFNm, PNEANet& Net, PNEANet& ReverseNet); //### build directed network graph
PNGraph LoadSIGraph(const TStr& InSIFNm, TIntIntVH& NIdInfTimeStepH, TIntIntVH& InfNIdInfoH);
void SaveSIGraph(const TStr& OutSIFNm, PNGraph InfCasc, TIntIntVH NIdInfTimeStepH, TIntIntVH InfNIdInfoH);
PNGraph RunSICascadeIC(const PNEANet& Net, const TIntV StartNodes, TIntIntVH& NIdInfTimeStepH, TIntIntVH& InfNIdInfoH, const TInt Ts); //### Independent Cascade model
PNGraph ReverseSIOrder(const PNEANet& G, const TIntH& ObservedNIdH, const TInt StartNode);
PNGraph ReverseSI(const PNEANet& G, const TInt StartNode);
TIntH Wavefront(const PNEANet& ReverseNet, const TIntH& ObservedNIdH, const TInt ObservedNId);
void getInfStat(const PNGraph& InfCasc, TIntH& ReverseNIdInfStat);
void MLCalc(TIntFltH& NIdStat, const TVec<TIntH>& TotalReverseStat);
TIntH TimeStepQ(const TIntIntVH& NIdInfTimeStepH, const TIntIntVH& InfNIdInfoH, const int TimeStep);
double reverseWeightVal(const TFlt weightVal);
double FltRnd();
int IntRnd();
void SaveQ(const TStr& OutSIFNm, const TIntH& SelectedNIdH, const TInt& InfSz);
void SaveCandidate(const TStr& OutSIFNm, const TIntFltH& NIdStat);
void SaveKPick(const TStr& OutSIFNm, const TIntV& K_Pick_Result, const TStr& Name);
void RandomKPick(const PNGraph& Sspct, TIntV& K_Pick_Result);
void TopKPickML(const TIntFltH& NIdStat, TIntV& K_Pick_Result);
void TopKPickInDeg(const PNGraph& Sspct, TIntV& K_Pick_Result);
void GreedyKPick(const PNGraph& Sspct, const TIntFltH& NIdStat, const TStr& flag, TIntV& K_Pick_Result);
TInt GreedyGetTopInDeg(const PNGraph& Sspct, TIntH& Covered);
TInt GreedyGetTopInDegML(const PNGraph& Sspct, const TIntFltH& NIdStat, TIntH& Covered);
PNGraph SuspectsCandidateGraph(const PNEANet& Net, const TIntFltH& NIdStat);
void NodeToHop(const PNEANet& Net, const TInt InfSrcNode, TIntV& K_Pick_Result);
void ApplyMinVal(TIntV& K_Pick_Result);
void ShiftMLToPositive(TIntFltH& NIdStat);
void PrintV(const TIntV& v, const TStr& s);
//
TFlt MinPenaltyPoint(const PNEANet& Net, const TIntFltH& NIdStat, const TIntV& KSet);
void PenaltyBasedKPick(const PNEANet& Net, const TIntFltH& NIdStat, TIntV& KSet);
void RdDummy (TIntFltH& NIdStat); 

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nTimeCascades. build: %s, %s. Start Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", INPUTFILE, "Input directed graph");
  if (DB==4) { printf("Loading %s\n", InFNm.CStr()); }// load input graph
  TFIn InFile(InFNm);
  PNEANet Net, ReverseNet;
  if (DB==4) { printf("Building Network...\n"); } 
  BuildNetwork(InFNm, Net, ReverseNet);
  TIntFltH NIdStat; RdDummy(NIdStat);
  NIdStat.SortByDat(false);
  printf("Output List: \n");
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    printf("<%d:%.4f>\n", (int)I.GetKey(), (double)I.GetDat());
  }
  printf("KSet: ");
  TIntV KSet; PenaltyBasedKPick(Net, NIdStat, KSet);
  for (int i=0; i<KSet.Len(); i++) {
    printf("%d, ", (int)KSet[i]);
  }  
  printf("\n");
  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}

//### build a directed network graph with weight value on each edge
void BuildNetwork(const TStr& InFNm, PNEANet& Net, PNEANet& ReverseNet) {
  TSsParser Ss(InFNm, ssfWhiteSep, true, true, true);
  Net.Clr(); ReverseNet.Clr();
  Net = TNEANet::New();
  ReverseNet = TNEANet::New();
  int SrcNId, DstNId, EId, ReverseEId; double WeightVal;
  while (Ss.Next()) {
    if (! Ss.GetInt(0, SrcNId) || ! Ss.GetInt(1, DstNId)) { continue; }
    Ss.GetFlt(2, WeightVal);
    if (! Net->IsNode(SrcNId)) { 
      Net->AddNode(SrcNId); //### add source node if it does not exist on the net
      ReverseNet->AddNode(SrcNId);
    }
    if (! Net->IsNode(DstNId)) {
      Net->AddNode(DstNId); //### add destination node if it does not exist on the net
      ReverseNet->AddNode(DstNId);
    }
    EId = Net->AddEdge(SrcNId, DstNId); //### connect the two nodes
    ReverseEId = ReverseNet->AddEdge(DstNId, SrcNId);
    Net->AddFltAttrDatE(EId, WeightVal, TSnap::CapAttrName); //### add weight value b/w nodes
    ReverseNet->AddFltAttrDatE(ReverseEId, WeightVal, TSnap::CapAttrName);
  }
  Net->Defrag(); //### defragment the net (compact the memory)
  ReverseNet->Defrag();
}

//### build a directed network graph with weight value on each edge
void BuildNetworkRandomWeight(const TStr& InFNm, PNEANet& Net, PNEANet& ReverseNet) {
  TSsParser Ss(InFNm, ssfWhiteSep, true, true, true);
  Net.Clr(); ReverseNet.Clr();
  Net = TNEANet::New();
  ReverseNet = TNEANet::New();
  int SrcNId, DstNId, EId, ReverseEId;
  double WeightVal;
  while (Ss.Next()) {
    if (! Ss.GetInt(0, SrcNId) || ! Ss.GetInt(1, DstNId)) { continue; }
    if (! Net->IsNode(SrcNId)) { 
      Net->AddNode(SrcNId); //### add source node if it does not exist on the net
      ReverseNet->AddNode(SrcNId);
    }
    if (! Net->IsNode(DstNId)) {
      Net->AddNode(DstNId); //### add destination node if it does not exist on the net
      ReverseNet->AddNode(DstNId);
    }
    EId = Net->AddEdge(SrcNId, DstNId); //### connect the two nodes
    ReverseEId = ReverseNet->AddEdge(DstNId, SrcNId);
    double const val[3] = {0.1,0.01,0.001}; //strong, medium, weak
    WeightVal = val[IntRnd()%3];
    Net->AddFltAttrDatE(EId, WeightVal, TSnap::CapAttrName); //### add weight value b/w nodes
    ReverseNet->AddFltAttrDatE(ReverseEId, WeightVal, TSnap::CapAttrName);
  }
  Net->Defrag(); //### defragment the net (compact the memory)
  ReverseNet->Defrag();
}
  
//###
PNGraph LoadSIGraph(const TStr& InSIFNm, TIntIntVH& NIdInfTimeStepH, TIntIntVH& InfNIdInfoH) {
  const TStr TSH = "TSH-" + InSIFNm;
  TSsParser SsTS(TSH, ssfWhiteSep, true, true, true);
  while (SsTS.Next()) {
    TIntV v;
    for (int i=1; i<SsTS.Len(); i++) {
      v.Add(atoi(SsTS[i]));
    }
    NIdInfTimeStepH.AddDat(atoi(SsTS[0]), v);
  }
  
  const TStr NIdInfo = "NIdInfo-" + InSIFNm;
  TSsParser SsInf(NIdInfo, ssfWhiteSep, true, true, true);
  while (SsInf.Next()) {
    TIntV v; int NId, TimeStep, InfSrc;
    if (SsInf.GetInt(0, NId) && SsInf.GetInt(1, TimeStep) && SsInf.GetInt(2, InfSrc)) {
      v.Add(TimeStep); v.Add(InfSrc);
      InfNIdInfoH.AddDat(NId, v);
    }
  }
  
  PNGraph Casc = TNGraph::New();
  TSsParser SsGraph(InSIFNm, ssfWhiteSep, true, true, true);
  int Src, Dst;
  while (SsGraph.Next()) {
    if (SsGraph.GetInt(0, Src) && SsGraph.GetInt(1, Dst)) { 
      if (!Casc->IsNode(Src)) {
        Casc->AddNode(Src);
      }
      if (!Casc->IsNode(Dst)) {
        Casc->AddNode(Dst);
      }
      Casc->AddEdge(Src, Dst);
    }
  }
  return Casc;
}

//###
void SaveSIGraph(const TStr& OutSIFNm, PNGraph InfCasc, TIntIntVH NIdInfTimeStepH, TIntIntVH InfNIdInfoH) {
  const TStr Inf = "si-" + OutSIFNm;
  TSnap::SaveEdgeList(InfCasc, Inf.CStr(), "Infected Directed Graph Using SI Model");

  const TStr TSH = "TSH-" + OutSIFNm;
  FILE* F = fopen(TSH.CStr(), "w");
  fprintf(F, "#Time Step\tInfNId\n");
  for (THash<TInt, TIntV>::TIter I = NIdInfTimeStepH.BegI(); I != NIdInfTimeStepH.EndI(); I++) {
    fprintf(F, "%d\t", (int)I.GetKey());
    TIntV v = I.GetDat();
    for (int i=0; i<v.Len(); i++) {
      fprintf(F, "%d ", (int)v[i]);
    }
    fprintf(F, "\n");
  }
  fclose(F);

  const TStr NIdInfo = "NIdInfo-" + OutSIFNm;
  F = fopen(NIdInfo.CStr(), "w");
  fprintf(F, "#NId\tTime Step\tInfSrc\n");
  for (THash<TInt, TIntV>::TIter I = InfNIdInfoH.BegI(); I != InfNIdInfoH.EndI(); I++) {
    fprintf(F, "%d\t", (int)I.GetKey());
    TIntV v = I.GetDat();
    for (int i=0; i<v.Len(); i++) {
      fprintf(F, "%d ", (int)v[i]);
    }
    fprintf(F, "\n");
  }
  fclose(F);
}

// Independent Cascade model (IC) 
// simulate SI model cascade using infection probability Threshold until stop 
PNGraph RunSICascadeIC(const PNEANet& G, const TIntV StartNodes, TIntIntVH& NIdInfTimeStepH, TIntIntVH& InfNIdInfoH, const TInt Ts) {
  if (DB==5) { printf("Running <RunSICascadeIC>\n"); }
  const double Epsilon = 0.0001; 
  PNGraph Casc = TNGraph::New();
  TIntH NIdActiveFlgH; //### flg: 0:inactive, 1:newly active, 2:active, etc.
  TIntQ Q;
  int timeStep = 0;
  TIntV InfNId, NIdInfo;
  for(int j=0; j<StartNodes.Len(); j++) { //### push all start nodes to queue
    const int StNId = StartNodes[j];
    Casc->AddNode(StNId);
    NIdActiveFlgH.AddDat(StNId, 0);
    Q.Push(StNId);
    InfNId.Add(StNId);
    NIdInfo.Add(timeStep); NIdInfo.Add(StNId);
    InfNIdInfoH.AddDat(StNId, NIdInfo);
    NIdInfo.Clr();
  }
  NIdInfTimeStepH.AddDat(timeStep, InfNId);
  while (!Q.Empty()) { 
    InfNId.Clr(); NIdInfo.Clr();
    if (timeStep==Ts) {
      return Casc;
    } else { timeStep += 1; }
    const TNEANet::TNodeI NI = G->GetNI(Q.Top()); Q.Pop();
    const int NodeId = NI.GetId();
    const int activeFlg = NIdActiveFlgH.GetDat(NodeId)+1;
    NIdActiveFlgH.AddDat(NodeId, activeFlg); //### set node active/infected
    int inActiveNeighbor = 0;
    for (int i = 0; i < NI.GetOutDeg(); i++) {
      NIdInfo.Clr();
      const int NeighborNId = NI.GetOutNId(i);
      if (Casc->IsNode(NeighborNId)) { continue; } //### node is already active/infected
      const int EId = G->GetEId(NodeId, NeighborNId);
      const double weightVal = G->GetFltAttrDatE(EId, TSnap::CapAttrName) * TMath::Power(0.75, activeFlg);
      if (weightVal < Epsilon) { continue; } //
      if (FltRnd() <= weightVal) {
        Casc->AddNode(NeighborNId); //### new active/infected node
        Casc->AddEdge(NodeId, NeighborNId);
        Q.Push(NeighborNId);
        NIdActiveFlgH.AddDat(NeighborNId, 0);
        InfNId.Add(NeighborNId);
        NIdInfo.Add(timeStep); NIdInfo.Add(NodeId);
        InfNIdInfoH.AddDat(NeighborNId, NIdInfo);
      } else {
        inActiveNeighbor += 1; //### count inactive neighbor here
      }
    }
    if (inActiveNeighbor){ //### check if there is still inactive neighbor
      Q.Push(NodeId);
    }
    if (InfNId.Len()!=0) { //### check if there is infected node
      NIdInfTimeStepH.AddDat(timeStep, InfNId);
    }
  }
  return Casc;
}

//### No infection order 
PNGraph ReverseSI(const PNEANet& G, const TInt StartNodes) {
  if(DB==5) { printf("in_ReverseSI\n"); }
  PNGraph Casc = TNGraph::New();
  Casc->AddNode(StartNodes);
  TIntQ Q; Q.Push(StartNodes);
  while (!Q.Empty()) {
    const TNEANet::TNodeI NI = G->GetNI(Q.Top()); Q.Pop();
    const int NodeId = NI.GetId();
    for (int i = 0; i < NI.GetOutDeg(); i++) {
      const int NeighborNId = NI.GetOutNId(i);
      if (Casc->IsNode(NeighborNId)) { continue; } //### node is already active/infected
      const int EId = G->GetEId(NodeId, NeighborNId); IAssert(EId>=0);
      const double weightVal = reverseWeightVal(G->GetFltAttrDatE(EId, TSnap::CapAttrName));
      if (FltRnd() <= weightVal) {
        Casc->AddNode(NeighborNId);
        Casc->AddEdge(NodeId, NeighborNId);
        Q.Push(NeighborNId);
      }
    }
  }
  if(DB==5) { printf("out_ReverseSI\n"); }
  return Casc;
}

//### Consider the order of infection
PNGraph ReverseSIOrder(const PNEANet& G, const TIntH& ObservedNIdH, const TInt StartNode) {
  if(DB==5) { printf("in_ReverseSI\n"); }
  TInt StartNodeTS = ObservedNIdH.GetDat(StartNode);
  TIntH BlacklistedNIdH, TSH; 
  PNGraph Casc = TNGraph::New();
  Casc->AddNode(StartNode);
  TIntQ Q; Q.Push(StartNode);
  while (!Q.Empty()) {
    const TNEANet::TNodeI NI = G->GetNI(Q.Top()); Q.Pop();
    const int NodeId = NI.GetId();
    TInt NodeIdTS = StartNodeTS;
    if(TSH.IsKey(NodeId)) { NodeIdTS = TSH.GetDat(NodeId); }
    for (int i = 0; i < NI.GetOutDeg(); i++) {
      const int NeighborNId = NI.GetOutNId(i);
      TInt NeighborNIdTS = NodeIdTS;
      if(ObservedNIdH.IsKey(NeighborNId)){
        if(ObservedNIdH.GetDat(NeighborNId)<NeighborNIdTS) {
          NeighborNIdTS = ObservedNIdH.GetDat(NeighborNId);
        } else {
          BlacklistedNIdH.AddDat(NeighborNId); continue;
        }
      } else if (TSH.IsKey(NeighborNId)) {
        if(TSH.GetDat(NeighborNId)<NeighborNIdTS) {
          NeighborNIdTS = TSH.GetDat(NeighborNId);
        } else {
          BlacklistedNIdH.AddDat(NeighborNId); continue;
        }
      }
      if (BlacklistedNIdH.IsKey(NeighborNId) || Casc->IsNode(NeighborNId)) { continue; } //### blacklisted or already active/infected
      const int EId = G->GetEId(NodeId, NeighborNId); IAssert(EId>=0);
      const double weightVal = reverseWeightVal(G->GetFltAttrDatE(EId, TSnap::CapAttrName));
      if ((FltRnd() <= weightVal)) {
        Casc->AddNode(NeighborNId); //### new active/infected node
        Casc->AddEdge(NodeId, NeighborNId);
        Q.Push(NeighborNId);
        TSH.AddDat(NeighborNId, NeighborNIdTS);// add/update timestep
      }
    }
  }
  if(DB==5) { printf("out_ReverseSI\n"); }
  return Casc;
}

//###
TIntH Wavefront(const PNEANet& ReverseNet, const TIntH& ObservedNIdH, const TInt ObservedNId) {
  if(DB==5) { printf("in_WF\n"); }
  TIntH ReverseNIdInfStat;
  TIntV StartNode; StartNode.Add(ObservedNId);
  TIntIntVH NIdInfTimeStepH, InfNIdInfoH; //### extra info of reverse spreading
  for (int simulationSize=0; simulationSize<SIMULATIONSZ; simulationSize++) {
#ifdef CONSIDER_ORDER
    PNGraph InfCasc = ReverseSIOrder(ReverseNet, ObservedNIdH, ObservedNId);
#else
    PNGraph InfCasc = ReverseSI(ReverseNet, ObservedNId);
#endif
    getInfStat(InfCasc, ReverseNIdInfStat); //### sum the occurence of nodes during SIMULATIONSZ run
  }
  printf("."); fflush(stdout);
  if(DB==5) { printf("in_WF\n"); }
  return ReverseNIdInfStat;
}

//### 
double reverseWeightVal(const TFlt weightVal) {
  double sum = 1.0;
  for (int i=0; i<20; i++) {
    sum = sum * (1 - (weightVal * TMath::Power(0.75, i)));
  }
  return(1-sum);
}

//### get stat of each reverse spreading node
void getInfStat(const PNGraph& InfCasc, TIntH& ReverseNIdInfStat) {
  if(DB==5) { printf("in_GetStat\n"); }
  for (TNGraph::TNodeI NI = InfCasc->BegNI(); NI < InfCasc->EndNI(); NI++) {
    if (ReverseNIdInfStat.IsKey(NI.GetId())) {
      ReverseNIdInfStat.AddDat(NI.GetId(), ReverseNIdInfStat.GetDat(NI.GetId())+1);
    } else {
      ReverseNIdInfStat.AddDat(NI.GetId(), 1);
    }
  }
  if(DB==5) { printf("out_GetStat\n"); }
}

//### calculate maximum likelihood 
void MLCalc(TIntFltH& NIdStat, const TVec<TIntH>& TotalReverseStat) {
  if(DB==5) { printf("in_ML\n"); }
  int index = 0;
  for (int i=1; i<TotalReverseStat.Len(); i++) {  //### find the smaller infected graph from reverse spreading
    if (TotalReverseStat[i].Len()<TotalReverseStat[index].Len()) { index = i; }
  }
  TIntH minInfCascH = TotalReverseStat[index];
  double probability = 0.0; int flg = 0;
  for (THash<TInt, TInt>::TIter I = minInfCascH.BegI(); I != minInfCascH.EndI(); I++) {
    for (int i=0; i<TotalReverseStat.Len(); i++) {
      if (TotalReverseStat[i].IsKey(I.GetKey())) {
        probability += TMath::Log(TotalReverseStat[i].GetDat(I.GetKey())) - TMath::Log(SIMULATIONSZ);
        flg = 1;
      } else {
        flg = 0; break;
      }
    }
    if(flg==1){ NIdStat.AddDat(I.GetKey(), probability); }
    probability = 0.0;
  }
  if(DB==5) { printf("out_ML\n"); }
}

//### Hash of (k,v) == (NId, InfTS)
TIntH TimeStepQ(const TIntIntVH& NIdInfTimeStepH, const TIntIntVH& InfNIdInfoH, const int TimeStep) {
  TIntV v; TIntQ Q; TIntH InfNId;
  for (THash<TInt, TIntV>::TIter I = NIdInfTimeStepH.BegI(); I != NIdInfTimeStepH.EndI(); I++) {
    if (I.GetKey()==0) { continue; }
    TIntV v = I.GetDat();
    for (int i=0; i<v.Len(); i++) {
      if (FltRnd()<Q_Prob) {
        if (Q.Len()>=Q_SZ) { Q.Pop(); }
        Q.Push(v[i]); 
      }
    }
  }
  for (int i=0; i<Q.Len(); i++) {
    InfNId.AddDat(Q[i], InfNIdInfoH.GetDat(Q[i])[0]);
  }
  return InfNId;
}

double FltRnd(){ //### 0<FltRnd()<1
  return (Rnd.GetUniDev());
}

int IntRnd(){ //### 
  return (Rnd.GetUniDevInt());
}

void SaveQ(const TStr& OutSIFNm, const TIntH& SelectedNIdH, const TInt& InfSz){
  const TStr QInfo = "Q-" + OutSIFNm;
  FILE* F = fopen(QInfo.CStr(), "a+");
  fprintf(F, "#NId\tTime Step\tInf-Graph: %d\n", (int)InfSz);
  for (THash<TInt, TInt>::TIter I = SelectedNIdH.BegI(); I != SelectedNIdH.EndI(); I++) {
    fprintf(F, "%d\t%d\n", (int)I.GetKey(), (int)I.GetDat());
  }
  fprintf(F, "\n");
  fclose(F); 
}

void SaveCandidate(const TStr& OutSIFNm, const TIntFltH& NIdStat){
  char fnm[128]; sprintf(fnm, "%lu-candidate-%s", (long)WallClockTime(), OutSIFNm.CStr());
  FILE* F = fopen(fnm, "a+");
  fprintf(F, "#NId\tML\n");
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    fprintf(F, "%d\t%.8f\n", (int)I.GetKey(), (double)I.GetDat());
  }
  fprintf(F, "\n");
  fclose(F); 
}

void SaveKPick(const TStr& OutSIFNm, const TIntV& K_Pick_Result, const TStr& Name){
  char fnm[128]; sprintf(fnm, "%s-%s", Name.CStr(), OutSIFNm.CStr());
  FILE* F = fopen(fnm, "a+");
  for (int i=0; i<K_Pick_Result.Len(); i++) {
    fprintf(F, "%d\t", (int)K_Pick_Result[i]);
  }
  fprintf(F, "\n");
  fclose(F); 
}

//reference to a vector of random pick nodes
//randomly pick k nodes from candidate set
void RandomKPick(const PNGraph& Sspct, TIntV& K_Pick_Result) {
  TRnd Random(IntRnd());
  while(K_Pick_Result.Len()<K_PICK) {
    const TInt NodeId = Sspct->GetRndNId(Random);
    if(K_Pick_Result.IsIn(NodeId)) { continue; }
    K_Pick_Result.Add(NodeId);
  }
} 

//select top k from candidate set using ML
void TopKPickML(const TIntFltH& NIdStat, TIntV& K_Pick_Result) {
  int count = 0;
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    K_Pick_Result.Add(I.GetKey());
    count = count + 1;
    if (count == K_PICK) { break; } 
  }
}

//select top k from candidate set using in-degree (not greedy)
void TopKPickInDeg(const PNGraph& Sspct, TIntV& K_Pick_Result) {
  TIntH DegStat; int count = 0;
  for (TNGraph::TNodeI NI = Sspct->BegNI(); NI < Sspct->EndNI(); NI++) {
    DegStat.AddDat(NI.GetId(), NI.GetInDeg());
  }
  DegStat.SortByDat(false); //biggest to smallest
  for(THash<TInt, TInt>::TIter tp = DegStat.BegI(); tp!=DegStat.EndI(); tp++) { //get biggest element (top)
    K_Pick_Result.Add(tp.GetKey());
    count = count + 1;
    if (count == K_PICK) { break; } 
  }
}

//general greedy pick with flag to switch inner func call
void GreedyKPick(const PNGraph& Sspct, const TIntFltH& NIdStat, const TStr& flag, TIntV& K_Pick_Result) {
  TIntH Covered; printf("Using <%s>\n", flag.CStr());
  for (int i=0; i<K_PICK; i++) {
    if (flag == "GreedyKPickInDeg") { //greedy pick but only consider top in-degree
      K_Pick_Result.Add(GreedyGetTopInDeg(Sspct, Covered));
    } else { //greedy pick but consider both in-degree and ML "GreedyKPickInDegML"
      K_Pick_Result.Add(GreedyGetTopInDegML(Sspct, NIdStat, Covered));
    }
  } 
}

TInt GreedyGetTopInDeg(const PNGraph& Sspct, TIntH& Covered) {
  TIntH DegStat; 
  for (TNGraph::TNodeI NI = Sspct->BegNI(); NI < Sspct->EndNI(); NI++) {
    int countInDeg = 0; const int NodeId = NI.GetId(); 
    if(Covered.IsKey(NodeId)) { continue; } 
    for (int i = 0; i < NI.GetInDeg(); i++) {
      if(!Covered.IsKey(NI.GetInNId(i))) { countInDeg = countInDeg + 1; } 
    }
    DegStat.AddDat(NodeId, countInDeg);
  }
  DegStat.SortByDat(false); //biggest to smallest
  THash<TInt, TInt>::TIter tp = DegStat.BegI(); //get biggest element (top)
  Covered.AddDat(tp.GetKey()); //mark node as convered
  return tp.GetKey();
}

TInt GreedyGetTopInDegML(const PNGraph& Sspct, const TIntFltH& NIdStat, TIntH& Covered) {
  TIntH Stat; 
  for (TNGraph::TNodeI NI = Sspct->BegNI(); NI < Sspct->EndNI(); NI++) {
    const int NodeId = NI.GetId(); double sumML = NIdStat.GetDat(NodeId);
    if(Covered.IsKey(NodeId)) { continue; } 
    for (int i = 0; i < NI.GetInDeg(); i++) {
      if(!Covered.IsKey(NI.GetInNId(i))) { 
        sumML = sumML + NIdStat.GetDat(NI.GetInNId(i));  
      } 
    }
    Stat.AddDat(NodeId, sumML);
  }
  Stat.SortByDat(false); //biggest to smallest
  THash<TInt, TInt>::TIter tp = Stat.BegI(); //get biggest element (top)
  Covered.AddDat(tp.GetKey()); //mark node as convered
  return tp.GetKey();
}

//convert hash table of candidate set into a directed graph 
PNGraph SuspectsCandidateGraph(const PNEANet& Net, const TIntFltH& NIdStat) {
  PNGraph Sspct = TNGraph::New();
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    const int NodeId = I.GetKey(); 
    if (!Sspct->IsNode(NodeId)) { Sspct->AddNode(NodeId); }
    const TNEANet::TNodeI NI = Net->GetNI(NodeId);
    for (int i = 0; i < NI.GetInDeg(); i++) {
      const int NeighborNId = NI.GetInNId(i);
      if(NIdStat.IsKey(NeighborNId)) { 
        if (!Sspct->IsNode(NeighborNId)) { Sspct->AddNode(NeighborNId); }
        Sspct->AddEdge(NeighborNId, NodeId);
      }
    }
    for (int i = 0; i < NI.GetOutDeg(); i++) {
      const int NeighborNId = NI.GetOutNId(i);
      if(NIdStat.IsKey(NeighborNId)) { 
        if (!Sspct->IsNode(NeighborNId)) { Sspct->AddNode(NeighborNId); }
        Sspct->AddEdge(NodeId, NeighborNId);
      }
    }
  }
  return Sspct;
}

void NodeToHop(const PNEANet& Net, const TInt InfSrcNode, TIntV& K_Pick_Result) {
  TBreathFS<PNEANet> BFS(Net);
  for (int i=0; i<K_Pick_Result.Len(); i++) {
    BFS.DoBfs(InfSrcNode, true, true, K_Pick_Result[i], TInt::Mx);
    K_Pick_Result[i] = BFS.GetHops(InfSrcNode, K_Pick_Result[i]);
  }
}

void ApplyMinVal(TIntV& K_Pick_Result) {
  TIntV tmp = K_Pick_Result; int min = tmp[0];
  for (int j=1; j<K_PICK; j++) {
    if(tmp[j]<min) {
      min = tmp[j];
    } else {
      K_Pick_Result[j] = min;
    }
  }
}

void ShiftMLToPositive(TIntFltH& NIdStat) {
  THash<TInt, TFlt>::TIter iter = NIdStat.EndI(); iter--;
  double minVal = -iter.GetDat() + 10.0;
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    const TInt NodeId = I.GetKey();
    NIdStat.AddDat(NodeId, NIdStat.GetDat(NodeId)+minVal);
  }
}

void PrintV(const TIntV& v, const TStr& s) {
  printf("%s", s.CStr());
  for (int i=0; i<(int)v.Len(); i++) {
    printf("%d ", (int)v[i]);
  } printf("\n");
}


void PenaltyBasedKPick(const PNEANet& Net, const TIntFltH& NIdStat, TIntV& KSet) { 
  TFlt minPen = TInt::Mx; TInt NId = -1;
  while (KSet.Len()<K_PICK) {
    for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
      if (KSet.IsIn(I.GetKey())) { continue; }
      TInt attemp = I.GetKey(); KSet.Add(attemp);
      TFlt Pen = MinPenaltyPoint(Net, NIdStat, KSet);
      KSet.DelLast();
      if (minPen>Pen) { minPen = Pen; NId = attemp; }
    }
    KSet.Add(NId);
  } 
}

TFlt MinPenaltyPoint(const PNEANet& Net, const TIntFltH& NIdStat, const TIntV& KSet) {
  TInt sp = TInt::Mx; TFlt sumPenalty = 0.0;
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    if (KSet.IsIn(I.GetKey())) { continue; }
    for (int i=0; i<KSet.Len(); i++) {
      TInt spath = TSnap::GetShortPath(Net, I.GetKey(), KSet[i]);
      if (sp>spath) { sp = spath; }
    }
    sumPenalty += sp * NIdStat.GetDat(I.GetKey());
    sp = TInt::Mx;
  }
  return sumPenalty;
}

void RdDummy (TIntFltH& NIdStat) {
  const TStr NIdInfo = "rdummy.txt";
  TSsParser SsInf(NIdInfo, ssfWhiteSep, true, true, true);
  while (SsInf.Next()) {
    int NId; double ml; 
    if (SsInf.GetInt(0, NId) && SsInf.GetFlt(1, ml)) {
      NIdStat.AddDat(NId, ml);
    }
  }
}
