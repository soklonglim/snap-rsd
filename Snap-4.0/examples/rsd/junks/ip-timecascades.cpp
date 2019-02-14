#include "stdafx.h"
#include <omp.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <dirent.h>

#define K_PICK 20
#define DBR 1
const TStr INPUTFILE = "Cit-HepTh.txt";

void NodeToHop(const PNEANet& Net, const TInt InfSrcNode, TIntV& K_Pick_Result);
void ApplyMinVal(TIntV& K_Pick_Result);
void ShiftMLToPositive(TIntFltH& NIdStat);
void PrintV(const TIntV& v, const TStr& s);
TFlt MinPenaltyPointPreProcessing(const TIntFltH& NIdStat, const TIntV& KSet, TStrH& pathMemory);
TFlt MinPenaltyPoint(const PNEANet& Net, const TIntFltH& NIdStat, const TIntV& KSet, TStrH& pathMemory);
void PenaltyBasedKPick(const TIntFltH& NIdStat, TIntV& KSet);
void MLToProbability(TIntFltH& NIdStat, TIntFltH& NIdStatProb);
void getSPListH(TStrH& pathMem);
void GetFiles(TStrV& fileList);

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nTimeCascades. build: %s, %s. Start Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  
  TStrV fileList; GetFiles(fileList);
  for (int i=0; i<fileList.Len(); i++) {
    printf("%3d: <%s>\n", i, fileList[i].CStr());
  }  
  
  
/*
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
*/
  TIntVV TotalKSet;
  TIntV penaltyBasedKPick;
  int NumThreads = 16;
  omp_set_num_threads(NumThreads);
  omp_lock_t lock;
  omp_init_lock(&lock);
  #pragma omp parallel private(penaltyBasedKPick)
  {
  #pragma omp for nowait
    for (int i=0; i<fileList.Len(); i++) {
      PenaltyBasedKPick(NIdStat, penaltyBasedKPick);
      NodeToHop(Net, InfSrcNode, penaltyBasedKPick);
      SaveKPick(INPUTFILE, penaltyBasedKPick, "regular-PenaltyBasedKPick");
      ApplyMinVal(penaltyBasedKPick);
      SaveKPick(INPUTFILE, penaltyBasedKPick, "min-PenaltyBasedKPick");
    }
    omp_set_lock(&lock);
    TotalKSet.Add(penaltyBasedKPIck);
    omp_unset_lock(&lock);
  }
  omp_destroy_lock(&lock);
/*  
  printf("PenaltyBasedKPick Step...\n");  
  TIntV penaltyBasedKPick; PenaltyBasedKPick(NIdStat, penaltyBasedKPick);
  NodeToHop(Net, InfSrcNode, penaltyBasedKPick);
  if (DBR==7) { PrintV(penaltyBasedKPick, "PenaltyBasedKPick: Regular Hop "); }
  SaveKPick(INPUTFILE, penaltyBasedKPick, "regular-PenaltyBasedKPick");
  ApplyMinVal(penaltyBasedKPick);
  if (DBR==7) { PrintV(penaltyBasedKPick, "Min Hop "); }
  SaveKPick(INPUTFILE, penaltyBasedKPick, "min-PenaltyBasedKPick");
*/
  return 0;
}

//###
void GetFiles(TStrV& fileList) {
  DIR *dir = opendir("."); struct dirent *ent;
  while ((ent = readdir(dir)) != NULL) {
    TStr file = ent->d_name;
    if ( file != "." && file != ".." ) {
      fileList.Add(file);
    }
  }
  closedir(dir);
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

void PenaltyBasedKPick(const TIntFltH& NIdStat, TIntV& KSet) {
  // printf("PenaltyBasedKPick\n");
  TFlt minPen = TInt::Mx; TInt NId = -1;
  TStrH pathMemory; getSPListH(pathMemory);
  
  while (KSet.Len()<K_PICK) {
    printf("KSet.Len: <%d>\n", (int)KSet.Len());
    for (THash<TInt, TFlt>::TIter I = NIdStat.BegI(); I != NIdStat.EndI(); I++) {
      if (KSet.IsIn(I.GetKey())) { continue; }
      KSet.Add(I.GetKey()); // attemp to add node to the set
      TFlt Pen = MinPenaltyPointPreProcessing(NIdStat, KSet, pathMemory);
      KSet.DelLast(); // remove the attemp node
      if (minPen>Pen) { minPen = Pen; NId = I.GetKey(); }
    }
    KSet.Add(NId);
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



