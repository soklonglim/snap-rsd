#include "stdafx.h"
#include <omp.h>
#include <time.h>
#include <math.h>

#define CONSIDER_ORDER

#define DB 51
#define DBR 1
//#define K_PICK 20
int K_PICK = 20;
#define Q_Prob 1
#define SIMULATIONSZ 1000

// idev --cpus-per-task=4 --partition=pddms --time=2:00:00
/* Cit-HepTh, CollegeMsg, email-Eu-core, twitter, Wiki-Vote, Email-Enron, soc-Epinions */
/*MxInDegNId* Cit-HepTh, CollegeMsg, email-Eu-core, twitter, Wiki-Vote, Email-Enron, soc-Epinions */
/*MxOutDegNId* Cit-HepTh<9905111>, CollegeMsg<9>, email-Eu-core<160>, twitter<115485051>, Wiki-Vote<2565>  CollegeMsg MidOutDeg<176>*/
const TStr INPUTFILE = "pen.txt";
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
void MLCalc(TIntFltH& NIdStat, const TVec<TIntH>& TotalReverseStat, TIntIntVH& NumInfect);
TIntH TimeStepQ(const TIntIntVH& NIdInfTimeStepH, const TIntIntVH& InfNIdInfoH, const int TimeStep);
double reverseWeightVal(const TFlt weightVal);
double FltRnd();
int IntRnd();
void SaveQ(const TStr& OutSIFNm, const TIntH& SelectedNIdH, const TInt& InfSz);
void SaveCandidate(const TStr& OutSIFNm, const TIntFltH& NIdStat, const TInt& InfSrcNode, const TIntIntVH& NumInfect);
void SaveKPick(const TStr& OutSIFNm, const TIntV& K_Pick_Result, const TStr& Name);
void RandomKPick(const PNGraph& Sspct, TIntV& K_Pick_Result);
void TopKPickML(const TIntFltH& NIdStat, TIntV& K_Pick_Result);
void TopKPickInDeg(const PNGraph& Sspct, TIntV& K_Pick_Result);
void TopKPickInDegEntireNetwork(const PNEANet& Net, const TIntFltH& NIdStat, TIntV& K_Pick_Result);
void GreedyKPick(const PNGraph& Sspct, const TIntFltH& NIdStat, const TStr& flag, TIntV& K_Pick_Result);
TInt GreedyGetTopInDeg(const PNGraph& Sspct, TIntH& Covered);
TInt GreedyGetTopInDegML(const PNGraph& Sspct, const TIntFltH& NIdStat, TIntH& Covered);
PNGraph SuspectsCandidateGraph(const PNEANet& Net, const TIntFltH& NIdStat);
void NodeToHop(const PNEANet& Net, const TInt InfSrcNode, TIntV& K_Pick_Result);
void ApplyMinVal(TIntV& K_Pick_Result);
void ShiftMLToPositive(TIntFltH& NIdStat);
void PrintV(const TIntV& v, const TStr& s);
TFlt MinPenaltyPointPreProcessing(const TIntFltH& NIdStat, const TIntV& KSet, TStrH& pathMemory);
TFlt MinPenaltyPoint(const PNEANet& Net, const TIntFltH& NIdStat, const TIntV& KSet, TStrH& pathMemory);
void PenaltyBasedKPick(const PNEANet& Net, const TIntFltH& NIdStat, TIntV& KSet);
void MLToProbability(TIntFltH& NIdStat, TIntFltH& NIdStatProb);
void getSPListH(TStrH& pathMem);
TInt DOBFS(const PNEANet& Net, const TInt& Src, const TInt& Dst);

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nTimeCascades. build: %s, %s. Start Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm; int InfSrcNode = -1;
  int index = 0; int srcIndex = -1; double ML = 0.0; int sus_sz = 0; TStr OutFNm;
  double t1, t2, t3, t4, t5;
  double t = WallClockTime(); double tt = t;
  Try
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", INPUTFILE, "Input directed graph");
  OutFNm = Env.GetIfArgPrefixStr("-o:", OUTPUTFILE, "Output directed graph"); 
  const TStr InSIFNm = Env.GetIfArgPrefixStr("-si:", INFECTEDFILE, "Input infected directed graph (SI Model)");
  const TInt TS = Env.GetIfArgPrefixInt("-ts:", 80, "Select Infected Node at Time Step");
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
    TRnd Random(IntRnd()); InfSrcNode = Net->GetRndNId(Random); 
    printf("InfSrcNode: %d, SIMU: %d\n", InfSrcNode, SIMULATIONSZ);
    //InfSrcNode = 5029; printf("InfSrcNode: %d, SIMU: %d\n", InfSrcNode, SIMULATIONSZ);
    TIntV StartNodes; StartNodes.Add(InfSrcNode);
    if (!Net->IsNode(InfSrcNode)) { printf("Node %d Does Not Exist!\n", InfSrcNode); return (0); }
    InfCasc = RunSICascadeIC(Net, StartNodes, NIdInfTimeStepH, InfNIdInfoH, TS);
    int inf_sz = InfCasc->GetNodes();
    //if (inf_sz<20) { printf("Infected Graph is too small...\n"); return (0); } 
    if (inf_sz<200) {
      Q_SZ = 20;
    } else {
      Q_SZ = InfCasc->GetNodes()/10 + 1;
    }
  }
  t2 = WallClockTime() - t; t = WallClockTime(); //load graph
  const TIntH SelectedNIdH = TimeStepQ(NIdInfTimeStepH, InfNIdInfoH, TS);
  TIntV SelectedNId;
  for (THash<TInt, TInt>::TIter I = SelectedNIdH.BegI(); I != SelectedNIdH.EndI(); I++) { SelectedNId.Add(I.GetKey()); }
  int NumThreads = SelectedNId.Len();
  if (NumThreads<1) { printf("No Infected Node...\n"); return (0); }
  if (DB) { 
    printf("Selected NId for ReverseSI: ");
    for (int i=0; i<SelectedNId.Len(); i++) { printf("%d, ", (int)SelectedNId[i]); } printf("\n");
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
/*  
  TVec<TIntH> TotalReverseStat;
  for(int i=0; i<NumThreads; i++) {
    TotalReverseStat.Add(Wavefront(ReverseNet, SelectedNIdH, SelectedNId[i]));
  }
*/
  TIntFltH NIdStat; TIntIntVH NumInfect;
  MLCalc(NIdStat, TotalReverseStat, NumInfect);
  t5 = WallClockTime() - t; //maximum likelihood calculation
  tt = WallClockTime() - tt; //total running time
  NIdStat.SortByDat(false);
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    if(index==-1) { index = 1; } else { index = index + 1; }
    if (I.GetKey() == InfSrcNode) { srcIndex = index; ML = (double)I.GetDat();}
  }
  printf("\nTotal Suspects: %d\n", NIdStat.Len()); sus_sz = NIdStat.Len();
  printf("Suspects Set is %1.2f%% of the entire network\n", ((double)NIdStat.Len()/Net->GetNodes())*100);
  printf("NId:%d\t%d\t%4.8f\n", InfSrcNode, srcIndex, ML);
  // if (sus_sz<K_PICK) { return(0); } 
  if (sus_sz<K_PICK) { K_PICK = sus_sz; printf("suspected size: %d\n", K_PICK); } 
  TBreathFS<PNEANet> BFS(Net); int cntr = 0; int hop_sz = sus_sz/100;
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    BFS.DoBfs(InfSrcNode, true, true, I.GetKey(), TInt::Mx);
    char hoplist1[128]; 
    sprintf(hoplist1, "echo \"%d\t->\t%d:\t%dhops\" >> hoplist-%s", 
                           InfSrcNode, (int)I.GetKey(), BFS.GetHops(InfSrcNode, I.GetKey()), InFNm.CStr()); 
    system(hoplist1);  if(++cntr==hop_sz) { break; }
  }
  SaveCandidate(InFNm, NIdStat, InfSrcNode, NumInfect);
  char saveTime[128]; sprintf(saveTime, "echo \"%4.8f\t%4.8f\t%4.8f\t%4.8f\t%4.8f\" >> tl-%s",
                                              t1, t2, t3, t4, t5, INPUTFILE.CStr()); system(saveTime);
  THash<TInt, TFlt>::TIter tp = NIdStat.BegI();
  char mlfile[128]; sprintf(mlfile, "echo \"%d\t%2.8f\t%d\t%2.8f\" >> ml-list-%s", 
                                 (int)tp.GetKey(), (double)tp.GetDat(), InfSrcNode, ML, INPUTFILE.CStr()); system(mlfile);
  SaveQ(InSIFNm, SelectedNIdH, InfCasc->GetNodes());
  
  printf("Generate Suspected Candidate Graph...\n");  
  const PNGraph Sspct = SuspectsCandidateGraph(Net, NIdStat);

  printf("RandomKPick Step...\n");  
  TIntV randomKPick; RandomKPick(Sspct, randomKPick); // random k pick
  NodeToHop(Net, InfSrcNode, randomKPick); 
  if (DBR==2) { PrintV(randomKPick, "RandomKPick: Regular Hop "); }
  SaveKPick(INPUTFILE, randomKPick, "regular-RandomKPick");
  ApplyMinVal(randomKPick); 
  if (DBR==2) { PrintV(randomKPick, "Min Hop "); }
  SaveKPick(INPUTFILE, randomKPick, "min-RandomKPick");

  printf("TopKPickML Step...\n");  
  TIntV topKPickML; TopKPickML(NIdStat, topKPickML); // top k ML
  NodeToHop(Net, InfSrcNode, topKPickML); 
  if (DBR==3) { PrintV(topKPickML, "TopKPickML: Regular Hop "); }
  SaveKPick(INPUTFILE, topKPickML, "regular-TopKPickML");
  ApplyMinVal(topKPickML);
  if (DBR==3) { PrintV(topKPickML, "Min Hop "); }
  SaveKPick(INPUTFILE, topKPickML, "min-TopKPickML");

  printf("TopKPickInDeg Step...\n");  
  TIntV topKPickInDeg; TopKPickInDeg(Sspct, topKPickInDeg); // top k in-degree consider only in suspected candidate set
  NodeToHop(Net, InfSrcNode, topKPickInDeg); 
  if (DBR==4) { PrintV(topKPickInDeg, "TopKPickInDeg: Regular Hop "); }
  SaveKPick(INPUTFILE, topKPickInDeg, "regular-TopKPickInDeg");
  ApplyMinVal(topKPickInDeg);
  if (DBR==4) { PrintV(topKPickInDeg, "Min Hop "); }
  SaveKPick(INPUTFILE, topKPickInDeg, "min-TopKPickInDeg");

  printf("TopKPickInDegEntireNetwork Step...\n");  
  TIntV topKPickInDegEntireNetwork; TopKPickInDegEntireNetwork(Net, NIdStat, topKPickInDegEntireNetwork); // top k in-degree of entire network
  NodeToHop(Net, InfSrcNode, topKPickInDegEntireNetwork); 
  if (DBR==4) { PrintV(topKPickInDegEntireNetwork, "TopKPickInDegEntireNetwork: Regular Hop "); }
  SaveKPick(INPUTFILE, topKPickInDegEntireNetwork, "regular-TopKPickInDegEntireNetwork");
  ApplyMinVal(topKPickInDegEntireNetwork);
  if (DBR==4) { PrintV(topKPickInDegEntireNetwork, "Min Hop "); }
  SaveKPick(INPUTFILE, topKPickInDegEntireNetwork, "min-TopKPickInDegEntireNetwork");

  printf("GreedyKPickInDeg Step...\n");  
  TIntV greedyKPickInDeg; GreedyKPick(Sspct, NIdStat, "GreedyKPickInDeg", greedyKPickInDeg); // greedy top k in-degree
  NodeToHop(Net, InfSrcNode, greedyKPickInDeg); 
  if (DBR==5) { PrintV(greedyKPickInDeg, "GreedyKPickInDeg: Regular Hop "); }
  SaveKPick(INPUTFILE, greedyKPickInDeg, "regular-GreedyKPickInDeg");
  ApplyMinVal(greedyKPickInDeg);
  if (DBR==5) { PrintV(greedyKPickInDeg, "Min Hop "); }
  SaveKPick(INPUTFILE, greedyKPickInDeg, "min-GreedyKPickInDeg");

  printf("GreedyKPickInDegML Step...\n");  
  ShiftMLToPositive(NIdStat);//convert ML to positive 
  TIntV greedyKPickInDegML; GreedyKPick(Sspct, NIdStat, "GreedyKPickInDegML", greedyKPickInDegML); // gree top k in-degree && ML
  NodeToHop(Net, InfSrcNode, greedyKPickInDegML); 
  if (DBR==6) { PrintV(greedyKPickInDegML, "GreedyKPickInDegML: Regular Hop "); }
  SaveKPick(INPUTFILE, greedyKPickInDegML, "regular-GreedyKPickInDegML");
  ApplyMinVal(greedyKPickInDegML);
  if (DBR==6) { PrintV(greedyKPickInDegML, "Min Hop "); }
  SaveKPick(INPUTFILE, greedyKPickInDegML, "min-GreedyKPickInDegML");

  printf("PenaltyBasedKPick Step...\n");  
  TIntV penaltyBasedKPick; PenaltyBasedKPick(Net, NIdStat, penaltyBasedKPick);
  NodeToHop(Net, InfSrcNode, penaltyBasedKPick);
  if (DBR==7) { PrintV(penaltyBasedKPick, "PenaltyBasedKPick: Regular Hop "); }
  SaveKPick(INPUTFILE, penaltyBasedKPick, "regular-PenaltyBasedKPick");
  ApplyMinVal(penaltyBasedKPick);
  if (DBR==7) { PrintV(penaltyBasedKPick, "Min Hop "); }
  SaveKPick(INPUTFILE, penaltyBasedKPick, "min-PenaltyBasedKPick");

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  char cmd[128]; sprintf(cmd, "echo \"%d\t%d\t%2.8f\t%d\t%8.4f\" >> %s", InfSrcNode, srcIndex, ML, sus_sz, tt, OutFNm.CStr()); system(cmd);
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
void MLCalc(TIntFltH& NIdStat, const TVec<TIntH>& TotalReverseStat, TIntIntVH& NumInfect) {
  if(DB==5) { printf("in_ML\n"); }
  int index = 0;
  for (int i=1; i<TotalReverseStat.Len(); i++) {  //### find the smaller infected graph from reverse spreading
    if (TotalReverseStat[i].Len()<TotalReverseStat[index].Len()) { index = i; }
  }
  TIntH minInfCascH = TotalReverseStat[index];
  TIntV first10Num; int cnt = 0;
  double probability = 0.0; int flg = 0;
  for (THash<TInt, TInt>::TIter I = minInfCascH.BegI(); I != minInfCascH.EndI(); I++) {
    for (int i=0; i<TotalReverseStat.Len(); i++) {
      if (TotalReverseStat[i].IsKey(I.GetKey())) {
        probability += TMath::Log(TotalReverseStat[i].GetDat(I.GetKey())) - TMath::Log(SIMULATIONSZ);
        if (cnt<10) { first10Num.Add(TotalReverseStat[i].GetDat(I.GetKey())); cnt = cnt + 1; }
        flg = 1;
      } else {
        flg = 0; break;
      }
    }
    if(flg==1){ NIdStat.AddDat(I.GetKey(), probability); NumInfect.AddDat(I.GetKey(), first10Num); }
    probability = 0.0; first10Num.Clr(); cnt = 0;
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

void SaveCandidate(const TStr& OutSIFNm, const TIntFltH& NIdStat, const TInt& InfSrcNode, const TIntIntVH& NumInfect){
  char fnm[128]; sprintf(fnm, "%lu-candidate-%s", (long)WallClockTime(), OutSIFNm.CStr());
  FILE* F = fopen(fnm, "a+");
  fprintf(F, "#NId\tML\n");
  fprintf(F, "%d\t%.8f\n", -101, (double)InfSrcNode);
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    fprintf(F, "%d\t%.8f\t", (int)I.GetKey(), (double)I.GetDat());
    TIntV tmp = NumInfect.GetDat(I.GetKey());
    fprintf(F, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 
    (int)tmp[0], (int)tmp[1], (int)tmp[2], (int)tmp[3], (int)tmp[4], (int)tmp[5], (int)tmp[6], (int)tmp[7], (int)tmp[8], (int)tmp[9]);
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

//select top k from candidate set using in-degree of entire network (not greedy)
void TopKPickInDegEntireNetwork(const PNEANet& Net, const TIntFltH& NIdStat, TIntV& K_Pick_Result) {
  TIntH DegStat; int count = 0;
  for (TNEANet::TNodeI NI = Net->BegNI(); NI < Net->EndNI(); NI++) {
    DegStat.AddDat(NI.GetId(), NI.GetInDeg());
  }
  DegStat.SortByDat(false); //biggest to smallest
  for(THash<TInt, TInt>::TIter tp = DegStat.BegI(); tp!=DegStat.EndI(); tp++) { //get biggest element (top)
    if (!NIdStat.IsKey(tp.GetKey())) { continue; }
    K_Pick_Result.Add(tp.GetKey());
    count = count + 1;
    if (count == K_PICK) { break; } 
  }
}

//general greedy pick with flag to switch inner func call
void GreedyKPick(const PNGraph& Sspct, const TIntFltH& NIdStat, const TStr& flag, TIntV& K_Pick_Result) {
  TIntH Covered; // printf("Using <%s>\n", flag.CStr());
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

void MyNodeToHop(const PNEANet& Net, const TInt InfSrcNode, TIntV& K_Pick_Result) {
  for (int i=0; i<K_Pick_Result.Len(); i++) {
    K_Pick_Result[i] = DOBFS(Net, InfSrcNode, K_Pick_Result[i]);
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
  double minVal = -iter.GetDat() + 1.0;
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
  // printf("PenaltyBasedKPick\n");
  TFlt minPen = TInt::Mx; TInt NId = -1;
  TStrH pathMemory; int flg = 0; const TStr exception = "Cit-HepTh.txt"; 
  if (INPUTFILE==exception) {
    printf("Loading Shortest Path List\n");
    printf("<%s> <%s>\n", INPUTFILE.CStr(), exception.CStr()); 
    getSPListH(pathMemory); flg = 1;
  }
  printf("preprocessing flg: %d\n", flg);
  while (KSet.Len()<K_PICK) {
    printf("KSet.Len: <%d>\n", (int)KSet.Len());
    for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
      if (KSet.IsIn(I.GetKey())) { continue; }
      KSet.Add(I.GetKey()); // attemp to add node to the set
      TFlt Pen;
      if (flg) {
        Pen = MinPenaltyPointPreProcessing(NIdStat, KSet, pathMemory);
      } else {
        Pen = MinPenaltyPoint(Net, NIdStat, KSet, pathMemory);
      } 
      KSet.DelLast(); // remove the attemp node
      if (minPen>Pen) { minPen = Pen; NId = I.GetKey(); }
    }
    KSet.Add(NId);
    printf("<NID: %d, Penalty: %d>\n", NId, minPen);
    minPen = TInt::Mx;
  }
}

TFlt MinPenaltyPointPreProcessing(const TIntFltH& NIdStat, const TIntV& KSet, TStrH& pathMemory) {
  // printf("MinPenaltyPointPreProcessing\n");
  TInt sp = TInt::Mx; TFlt sumPenalty = 0.0;
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    if (KSet.IsIn(I.GetKey())) { continue; }
    for (int i=0; i<KSet.Len(); i++) {
      int src = I.GetKey(); int dst = KSet[i]; int spath; char e[128]; 
      if (src<dst) { 
        sprintf(e, "%d-%d", src, dst); 
      } else { 
        sprintf(e, "%d-%d", dst, src); 
      }
      if (pathMemory.IsKey(e)) {
        spath = pathMemory.GetDat(e);
      } else { // no shortest path found
        spath = TInt::Mx; printf("No Shortest Path Found <%d<->%d>\n", src, dst);
      }
      if (sp>spath) { sp = spath; if (sp==1) { break; } }
    }
    sumPenalty += sp * NIdStat.GetDat(I.GetKey());
    sp = TInt::Mx;
  }
  return sumPenalty;
}

void MLToProbability(TIntFltH& NIdStat, TIntFltH& NIdStatProb) {
  // printf("MLToProbability\n");
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    NIdStatProb.AddDat(I.GetKey(), TMath::Power(10.0, I.GetDat()));
  }
}

void getSPListH(TStrH& pathMem) {
  const TStr SPFILE = "SPList" + INPUTFILE;
  TSsParser Ss(SPFILE, ssfWhiteSep, true, true, true);
  while(Ss.Next()) {
    int p; Ss.GetInt(1, p);
    const TStr e = Ss[0];
    pathMem.AddDat(e, p);
  }
}

TFlt MinPenaltyPoint(const PNEANet& Net, const TIntFltH& NIdStat, const TIntV& KSet, TStrH& pathMem) {
  // printf("MinPenaltyPoint\n");
  TInt sp = TInt::Mx; TFlt sumPenalty = 0.0;
  for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
    if (KSet.IsIn(I.GetKey())) { continue; }
    for (int i=0; i<KSet.Len(); i++) {
      int src = I.GetKey(); int dst = KSet[i]; int spath; char e[128]; 
      if (src<dst) { 
        sprintf(e, "%d-%d", src, dst); 
      } else { 
        sprintf(e, "%d-%d", dst, src); 
      }
      if (pathMem.IsKey(e)) {
        spath = pathMem.GetDat(e);
      } else { 
        spath = TSnap::GetShortPath(Net, src, dst);
        pathMem.AddDat(e, spath); 
      }
      if (sp>spath) { sp = spath; if (sp==1) { break; } }
    }
    sumPenalty += sp * NIdStat.GetDat(I.GetKey());
    sp = TInt::Mx;
  }
  return sumPenalty;
}

TInt DOBFS(const PNEANet& Net, const TInt& Src, const TInt& Dst) {
  if(Src==Dst) { return 0; }
  TIntQ Q1, Q2; // queue for level1 & level2 search
  Q1.Push(Src);
  while(!Q1.Empty()) { //leve1 search
    TNEANet::TNodeI NI = Net->GetNI(Q1.Top()); Q1.Pop();
    for (int i=0; i<NI.GetOutDeg(); i++) { //out-links
      if(NI.GetOutNId(i)==Dst) { return 1; }
      Q2.Push(NI.GetOutNId(i));
    }
    for (int i=0; i<NI.GetInDeg(); i++) { //in-links
      if(NI.GetInNId(i)==Dst) { return 1; }
      Q2.Push(NI.GetInNId(i));
    }
  } //done level1 search

  while(!Q2.Empty()) { //level2 search
    TNEANet::TNodeI NI = Net->GetNI(Q2.Top()); Q2.Pop();
    for (int i=0; i<NI.GetOutDeg(); i++) { //out-links
      if(NI.GetOutNId(i)==Dst) { return 2; }
    }
    for (int i=0; i<NI.GetInDeg(); i++) { //in-links
      if(NI.GetInNId(i)==Dst) { return 2; }
    }
  } //done level2 search
  return 3;
}


