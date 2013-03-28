#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>

#define GZ_BUF_SIZE 1048576

using namespace std;

struct alignment {
	int total_len;
	int total_mis;
	int total_gap;
	int read_start;
	int read_end;
	int adpt_start;
	int adpt_end;
	float mis_rate;
};

int trace_back (int *DPscore, int *direction, vector<char>& aligni, vector<char>& alignj);
void split (string &strLine, vector<string>& tokens, const char* delim) ;
int max_score (int subscore, int gapiscore, int gapjscore, int &maxscore, int &direction);
alignment local_alignment ( string *read, string *adpt,map< char, map<char,int> > scoreMatrix,int gapPenalty);
