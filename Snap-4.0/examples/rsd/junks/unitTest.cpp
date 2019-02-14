#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <thread>
#include <ctime>
#include <chrono>

#define T_SZ 25
using namespace std;

// void mySysCall(int tid, int TimeStep) {
// 	// string id = std::to_string(tid);
// 	// string ts = std::to_string(TimeStep);
// 	// string cmd("./timecascades -ts:" + ts +" >> " + id + "bigoutput.txt\n");
// 	// system(cmd.c_str());
// 	// printf("%s", cmd.c_str());
// 	system("./timecascades -si:bigtest1.txt -ts:80");
// }

void mySysCall(int tid, int TimeStep) {
	system("./timecascades-Wiki-Vote");
}

int main(int argc, char** argv) {
	time_t now = time(0);
	char* dt = ctime(&now);


	int batchNum = 0;
	if (argc>1) {
		batchNum = atoi(argv[1]);
	} 

	thread t[T_SZ];
	for (int i=0; i<T_SZ; i++) {
		t[i] = thread(mySysCall, batchNum, 80);
		this_thread::sleep_for (std::chrono::seconds(1));
	}
	
	for (int i=0; i<T_SZ; i++) {
		t[i].join();
	}

	cout << "\nStart: " << dt << endl;
	now = time(0);
	dt = ctime(&now);
	cout << "Stop: " << dt << endl;
	return 0;
}


