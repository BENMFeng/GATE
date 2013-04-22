
#include "local_alignment.h"

using namespace std;


int trace_back (int *DPscore, int *direction,vector<char>& aligni,vector<char>& alignj);
int max_score (int subscore, int gapiscore, int gapjscore, int &maxscore, int &direction);

string seqi,seqj;
int ilength,jlength;
int posi=1;
int posj=1;
int mismatch_base;
int gap_length;
int total_len;



alignment local_alignment ( string *read, string *adpt,map< char, map<char,int> > scoreMatrix,int gapPenalty) {

	seqi = *read;
	seqj = *adpt;
	ilength = seqi.length();
	jlength = seqj.length();
	alignment alignment;
	mismatch_base=0;
	gap_length=0;
	total_len=0;
	/*		construct DP matrix		*/
	int *DPscore=new int[(ilength+1)*(jlength+1)+1];
	int *direction=new int[(ilength+1)*(jlength+1)+1];
	DPscore[0]=0;
	direction[0]=0;
	DPscore[1]=0;
	direction[1]=0;
	vector<char> aligni,alignj;
	for (int j=1; j<jlength+1; j++) {
		DPscore[j*(ilength+1)+(0+1)]=0;
		direction[j*(ilength+1)+(0+1)]=1;
		//0:from (i-1,j-1) to (i,j); 1:from (i-1,j) to (i,j); 2:from (i,j-1) to (i,j);
	}

	for (int i=1; i<ilength+1; i++) {
		DPscore[0*(ilength+1)+(i+1)]=0;
		direction[0*(ilength+1)+(i+1)]=2;
	}

	for (int i=1; i<ilength+1; i++) {
		for (int j=1; j<jlength+1; j++) {
			int maxi=0,maxj=0;

			int subsScore=DPscore[(j-1)*(ilength+1)+(i-1)+1]+scoreMatrix[seqi[i-1]][seqj[j-1]];
			//(i-1,j-1)->(i,j)

			for(int k=0;k<=j-1;k++)
				if (maxi<DPscore[k*(ilength+1)+(i+1)])
					maxi=DPscore[k*(ilength+1)+(i+1)];
			int gapiScore=maxi+gapPenalty;//(i,j-1)->(i,j)

			for(int k=0;k<=i-1;k++)
				if (maxj<DPscore[j*(ilength+1)+(k+1)])
					maxj=DPscore[j*(ilength+1)+(k+1)];
			int gapjScore=maxj+gapPenalty;//(i-1,j)->(i,j)
			max_score(subsScore,gapiScore,gapjScore,DPscore[j*(ilength+1)+(i+1)],direction[j*(ilength+1)+(i+1)]);

			if (DPscore[0]<=DPscore[j*(ilength+1)+(i+1)]) {
				posi=i;posj=j;//(posi,posj):position of the max score in the DP matrix;
				DPscore[0]=DPscore[j*(ilength+1)+(i+1)];//max score
			}
		}
	}
	alignment.read_end = posi-1;
	alignment.adpt_end = posj-1;

/*******	find the best alignment according to the DP matrix	*********/

	trace_back(DPscore,direction,aligni,alignj);

	alignment.read_start = posi;
	alignment.adpt_start = posj;
	alignment.total_len = total_len;
	alignment.total_mis = mismatch_base;
	alignment.total_gap = gap_length;

	reverse(aligni.begin(),aligni.end());
	reverse(alignj.begin(),alignj.end());

// 	/*		print the alignment		*/
// 	cout<<total_len<< " "<<mismatch_base<<" "<<gap_length<<endl;
// 	vector<char>::iterator myiter;
// 	myiter = aligni.begin();
// 	cout<<" "<<"seq1:";
// 	for ( int i = 0; myiter != aligni.end(); ++myiter, ++i )
// 	{
// 		cout << *myiter;
// 	}
// 	cout << endl;
// 	myiter = alignj.begin();
// 	cout<<" "<<"seq2:";
// 	for ( int i = 0; myiter != alignj.end(); ++myiter, ++i )
// 	{
// 		cout << *myiter;
// 	}
// 	cout << endl;

        delete DPscore;
        delete direction;

	return (alignment);



}

//dealing global varibles
int trace_back (int *DPscore, int *direction,vector<char>& aligni,vector<char>& alignj)
{
	if (DPscore[posj*(ilength+1)+(posi+1)]==0)
	{
		return(0);//stop,when meet the 0;
	}else{
		total_len++;
		switch (direction[posj*(ilength+1)+(posi+1)]) {
			case 0:{
				if(DPscore[posj*(ilength+1)+(posi+1)]<DPscore[(posj-1)*(ilength+1)+posi]) {
					mismatch_base++;
				}
				aligni.push_back(seqi[posi-1]);
				alignj.push_back(seqj[posj-1]);
				posi--;
				posj--;
				break;
			}
			case 1:aligni.push_back('-'); alignj.push_back(seqj[posj-1]); posj--;gap_length++; break;
			case 2:aligni.push_back(seqi[posi-1]); alignj.push_back('-'); posi--; gap_length++;break;
		}
		if (posi > 0 && posj > 0){
			trace_back(DPscore,direction,aligni,alignj);
		}
	}

}

int max_score (int subscore, int gapiscore, int gapjscore, int &maxscore, int &direction) {
	if (subscore<0)
		subscore=0;//if the score is negative,give it the value of 0;
	if (gapiscore<0)
		gapiscore=0;
	if (gapjscore<0)
		gapjscore=0;
	if (subscore >= gapiscore)
	{
		maxscore = subscore;
		direction = 0;
	}else{
		maxscore = gapiscore;
		direction = 1;
	}

	if (maxscore < gapjscore)
	{
		maxscore = gapjscore;
		direction = 2;
	}

}

void split(string &strLine, vector<string>& tokens, const char* delim)
{
	int count = 0;
	for(;;) {
		//erase delimiter
	   int i = strLine.find_first_not_of(delim);
	   if(i == -1)
	    	break;
	   strLine.erase(0, i);

	   i = strLine.find_first_of(delim);
	   if(i == -1) {
		    tokens.push_back(strLine);
		    break;
	   } else {
		    string token = strLine.substr(0, i);
		    strLine.erase(0, i);
		    tokens.push_back(token);
	   }
	}
}



