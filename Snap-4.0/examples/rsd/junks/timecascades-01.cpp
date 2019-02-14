#include "stdafx.h"
#include <omp.h>
#include <time.h>
#include <math.h>

#define CONSIDER_ORDER

#define DB 51
// #define Q_SZ 40
#define Q_Prob 0.15
#define SIMULATIONSZ 10000
// #define InfSrcNode 23 // source of infected node 

// idev --cpus-per-task=4 --partition=pddms --time=2:00:00
/* Cit-HepTh, CollegeMsg, email-Eu-core, Wiki-Vote, Email-Enron, */
const TStr INPUTFILE = "soc-pokec-relationships.txt";
const TStr OUTPUTFILE = "output-" + INPUTFILE;
const TStr INFECTEDFILE = "si-" + INPUTFILE;

int Q_SZ = 0;

TRnd Rnd(0);
void BuildNetwork(const TStr& InFNm, PNEANet& Net, PNEANet& ReverseNet); //### build directed network graph
PNGraph LoadSIGraph(const TStr& InSIFNm, TIntIntVH& NIdInfTimeStepH, TIntIntVH& InfNIdInfoH);
void SaveSIGraph(const TStr& OutSIFNm, PNGraph InfCasc, TIntIntVH NIdInfTimeStepH, TIntIntVH InfNIdInfoH);
PNGraph RunSICascadeIC(const PNEANet& Net, const TIntV StartNodes, TIntIntVH& NIdInfTimeStepH, TIntIntVH& InfNIdInfoH, const TInt Ts); //### Independent Cascade model
PNGraph ReverseSIOrder(const PNEANet& G, const TIntH& ObservedNIdH, const TInt StartNode);
PNGraph ReverseSI(const PNEANet& G, const TInt StartNode);
TIntH Wavefront(const PNEANet& ReverseNet, const TIntH& ObservedNIdH, const TInt ObservedNId);
void getInfStat(const PNGraph& InfCasc, TIntH& ReverseNIdInfStat);
void MLCalc(TIntFltH& NIdStat, const TVec<TIntH>& TotalReverseStat);
TIntH TimeStepQ(const TIntIntVH& NIdInfTimeStepH, const TIntIntVH& InfNIdInfoH, const int TimeStep);
float reverseWeightVal(const TFlt weightVal);
float FltRnd();
int IntRnd();
void SaveQ(const TStr& OutSIFNm, const TIntH& SelectedNIdH);

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nTimeCascades. build: %s, %s. Start Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm; int InfSrcNode = -1;
  int index = 0; int srcIndex = -1; float ML = 0.0; int sus_sz = 0; TStr OutFNm;
  double t1, t2, t3, t4, t5;
  double t = WallClockTime(); double tt = t;
  Try
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", INPUTFILE, "Input directed graph");
  OutFNm = Env.GetIfArgPrefixStr("-o:", OUTPUTFILE, "Output directed graph"); 
  const TStr InSIFNm = Env.GetIfArgPrefixStr("-si:", INFECTEDFILE, "Input infected directed graph (SI Model)");
  const TInt TS = Env.GetIfArgPrefixInt("-ts:", 200, "Select Infected Node at Time Step");
  if (DB==4) { printf("Loading %s\n", InFNm.CStr()); }// load input graph
  TFIn InFile(InFNm);
  PNEANet Net, ReverseNet;
  TIntIntVH NIdInfTimeStepH; //### (key, value) = (timestep, [list of all infected node at timestep])
  TIntIntVH InfNIdInfoH; //### (key, value) = (NId, [timestep, infSrcNId])
  if (DB==4) { printf("Building Network...\n"); } 
  BuildNetwork(InFNm, Net, ReverseNet);
  t1 = WallClockTime() - t; t = WallClockTime(); //build graph
  if (DB==4) { 
    printf("Printing Network to %s ...\n", OutFNm.CStr());
    TSnap::SaveEdgeListNet(Net, OutFNm.CStr(), NULL);
  }
  if (DB==4) { printf("Total Nodes:%d  Edges:%d\n", Net->GetNodes(), Net->GetEdges()); }
  PNGraph InfCasc; // Simulate Infected model
  if (InSIFNm != INFECTEDFILE) { // load graph from files
    InfCasc = LoadSIGraph(InSIFNm, NIdInfTimeStepH, InfNIdInfoH);
  } else { //simulate infected graph using SI model & save graph to files
    TRnd Random(IntRnd());
    InfSrcNode = Net->GetRndNId(Random); InfSrcNode = 28230;
    TIntV StartNodes; StartNodes.Add(InfSrcNode);
    if (!Net->IsNode(InfSrcNode)) { printf("Node %d Does Not Exist!\n", InfSrcNode); return (0); }
    InfCasc = RunSICascadeIC(Net, StartNodes, NIdInfTimeStepH, InfNIdInfoH, TS);
    // SaveSIGraph(InSIFNm, InfCasc, NIdInfTimeStepH, InfNIdInfoH);
    Q_SZ = InfCasc->GetNodes()/10 + 1;
  }
  printf("InfSrcNode: %d, InfGraphSize: %d\n", InfSrcNode, InfCasc->GetNodes());
  t2 = WallClockTime() - t; t = WallClockTime(); //load graph
  const TIntH SelectedNIdH = TimeStepQ(NIdInfTimeStepH, InfNIdInfoH, TS);
  TIntV SelectedNId;
  for (THash<TInt, TInt>::TIter I = SelectedNIdH.BegI(); I != SelectedNIdH.EndI(); I++) { SelectedNId.Add(I.GetKey()); }
  int NumThreads = SelectedNId.Len();
  if (NumThreads<1) { printf("No Infected Node...\n"); return (0); }
  if (DB) { 
    printf("Selected NId for ReverseSI: ");
    for (int i=0; i<SelectedNId.Len(); i++) { printf("%d, ", (int)SelectedNId[i]); } printf("\n"); //printf("\n<NId, TS>: ");
    // for (THash<TInt, TInt>::TIter I = SelectedNIdH.BegI(); I != SelectedNIdH.EndI(); I++) { printf("<%d,%d> ", (int)I.GetKey(), (int)I.GetDat()); } printf("\n");
  }
  t3 = WallClockTime() - t; t = WallClockTime(); //sampling selection
  for(int i=0; i<Q_SZ; i++) { printf("-"); } printf("\n");
  TVec<TIntH> TotalReverseStat;
  TIntH SingleInfH;
  omp_set_num_threads(NumThreads);
  omp_lock_t lock;
  omp_init_lock(&lock);
  #pragma omp parallel private(SingleInfH)
  {
  #pragma omp for nowait
    for(int i=0; i<NumThreads; i++){
      SingleInfH = Wavefront(ReverseNet, SelectedNIdH, SelectedNId[i]);
    }
    omp_set_lock(&lock);
    TotalReverseStat.Add(SingleInfH);
    omp_unset_lock(&lock);
  }
  omp_destroy_lock(&lock);
  t4 = WallClockTime() - t; t = WallClockTime(); //simulation
  TIntFltH NIdStat; 
  MLCalc(NIdStat, TotalReverseStat);
  t5 = WallClockTime() - t; //maximum likelihood calculation
  tt = WallClockTime() - tt; //total running time
  NIdStat.SortByDat(false);
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    index = index + 1;
    if (I.GetKey() == InfSrcNode) { srcIndex = index - 1; ML = (float)I.GetDat();}
  }
  printf("\nTotal Suspects: %d\n", NIdStat.Len()); sus_sz = NIdStat.Len();
  printf("Suspects Set is %1.2f%% of the entire network\n", ((float)NIdStat.Len()/Net->GetNodes())*100);
  printf("NId:%d\t%d\t%4.8f\n", InfSrcNode, srcIndex, ML);
  TBreathFS<PNEANet> BFS(Net); int cntr = 0; int hop_sz = sus_sz/100;
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    BFS.DoBfs(InfSrcNode, true, true, I.GetKey(), TInt::Mx);
    char hoplist1[128]; 
    sprintf(hoplist1, "echo \"%d\t->\t%d:\t%dhops\" >> hoplist-%s", 
                           InfSrcNode, (int)I.GetKey(), BFS.GetHops(InfSrcNode, I.GetKey()), InFNm.CStr());
    if(system(hoplist1));  if(++cntr==hop_sz) { break; }
  }
  char saveTime[128]; sprintf(saveTime, "echo \"%4.8f\t%4.8f\t%4.8f\t%4.8f\t%4.8f\" >> tl-%s",
                                              t1, t2, t3, t4, t5, INPUTFILE.CStr()); if(system(saveTime));
  THash<TInt, TFlt>::TIter tp = NIdStat.BegI();
  char mlfile[128]; sprintf(mlfile, "echo \"%d\t%2.8f\t%d\t%2.8f\" >> ml-list-%s", 
                                 (int)tp.GetKey(), (float)tp.GetDat(), InfSrcNode, ML, INPUTFILE.CStr()); if(system(mlfile));
  SaveQ(InSIFNm, SelectedNIdH);
  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  char cmd[128]; sprintf(cmd, "echo \"%d\t%d\t%2.8f\t%d\t%8.4f\" >> %s", InfSrcNode, srcIndex, ML, sus_sz, tt, OutFNm.CStr()); if(system(cmd));
  return 0;
}

//### build a directed network graph with weight value on each edge
void BuildNetwork(const TStr& InFNm, PNEANet& Net, PNEANet& ReverseNet) {
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
      const float weightVal = G->GetFltAttrDatE(EId, TSnap::CapAttrName) * TMath::Power(0.75, activeFlg);
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
      const float weightVal = reverseWeightVal(G->GetFltAttrDatE(EId, TSnap::CapAttrName));
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
      const float weightVal = reverseWeightVal(G->GetFltAttrDatE(EId, TSnap::CapAttrName));
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
float reverseWeightVal(const TFlt weightVal) {
  float sum = 1.0;
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
  float probability = 0.0; int flg = 0;
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

float FltRnd(){ //### 0<FltRnd()<1
  return (Rnd.GetUniDev());
}

int IntRnd(){ //### 
  return (Rnd.GetUniDevInt());
}

void SaveQ(const TStr& OutSIFNm, const TIntH& SelectedNIdH){
  const TStr QInfo = "Q-" + OutSIFNm;
  FILE* F = fopen(QInfo.CStr(), "a+");
  fprintf(F, "#NId\tTime Step\n");
  for (THash<TInt, TInt>::TIter I = SelectedNIdH.BegI(); I != SelectedNIdH.EndI(); I++) {
    fprintf(F, "%d\t%d\n", (int)I.GetKey(), (int)I.GetDat());
  }
  fprintf(F, "\n");
  fclose(F); 
}

