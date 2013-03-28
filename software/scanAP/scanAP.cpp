/*
 * scanAP.cpp

 *
 *  Created on: Oct 31, 2012
 *      Author: BENM binxiaofeng@gmail.com
 */
#include "local_alignment.h"
#include <iomanip>
using namespace std;

alignment get_alignment (string *read, int read_pos, string *adpt, int adapt_pos, int seed_len, int max_mis);
void usage();
void get_scoreMatrix (int &gapPenalty, map< char, map<char,int> > &scoreMatrix);

int main (int argc, char *argv[ ]) {

	int seed_len = 10;  //length of seed
	int max_mis = 3;    //the max mismatches
	int min_match = 10; //the min match length
	int opt;
	int ifrag = 0;
	int afrag = 0;
	int lfrag = 0;
	string readfile, adptfile; //input file
	string listfile, statfile; //output file
	alignment alignment;
	int gapPenalty;
	map< char, map<char,int> > scoreMatrix;
	while ((opt=getopt(argc, argv, "i:a:s:d:k:e:m:lh"))!=-1) {
		switch (opt) {
			case 'i': ifrag=1; readfile=optarg; break; //input file
			case 'a': afrag=1; adptfile=optarg; break; //input file
			case 's': statfile=optarg; break;          //output file
			case 'd': listfile=optarg; break;          //output file
			case 'k': seed_len=atoi(optarg); break;
			case 'e': max_mis=atoi(optarg); break;
			case 'm': min_match=atoi(optarg); break;
			case 'l': lfrag=1; break;
			case 'h': usage(); break;
			default:  usage();
		}
	}
	if (!(ifrag&&afrag)) {
		usage();
	}
	if (lfrag) {
		get_scoreMatrix (gapPenalty, scoreMatrix);
	}

	fstream fastq_file, adpt_file; //input file
	fastq_file.open (readfile.c_str(), ios::in);
	adpt_file.open (adptfile.c_str(), ios::in);
	if (fastq_file==NULL || adpt_file==NULL) {
		cout<<"Cann't open files!"<<endl;
		usage();
	}

	fstream list_file, stat_file; //output file
	list_file.open (listfile.c_str(), ios::out);
	stat_file.open (statfile.c_str(), ios::out);
	list_file<<"#reads_id   reads_len   reads_start   reads_end   adapter_id   adapter_len   adapter_start   adapter_end   align_len   mismatch   gap"<<endl;
	stat_file<<"adapter_id\tpolluted_reads\tempty_reads\tadapter_sequence"<<endl;

	string read_id, read, qual_id, qual;
	string adpt_id, adpt;
	int read_pos; //start position in reads
	int adpt_pos; //start position in adapter
	int read_num = 0;
	int adpt_num = 0;
	int empty_num = 0;

	while (getline (adpt_file, adpt_id, '\n')) {
		adpt_id.erase(0,1);
		getline (adpt_file, adpt, '\n');
		int adpt_len = adpt.length();

		read_num=0;
		adpt_num=0;
		empty_num = 0;

		while (getline (fastq_file, read_id, '\n')) {
			read_id.erase(0,1);
			getline (fastq_file, read, '\n');
			int read_len = read.length();
			getline (fastq_file, qual_id, '\n');
			getline (fastq_file, qual, '\n');
			read_num++;

			int find;
			for (adpt_pos=0; adpt_pos<=adpt_len-seed_len; adpt_pos++) {
				for (read_pos=0; read_pos<=read_len-seed_len; read_pos++) {
					find = 1;
					for (int i=0; i<seed_len; i++) {
						if (adpt[adpt_pos+i]!=read[read_pos+i]) {
							find = 0;
							break;
						}
					}
					if (find) {
						break;
					}
				}
				if (find) {
					break;
				}
			}

			if (find) {
				if (!lfrag) {
					alignment = get_alignment (&read, read_pos, &adpt, adpt_pos, seed_len, max_mis);
				}
				else{
					alignment = local_alignment (&read, &adpt, scoreMatrix, gapPenalty);
				}
				if (alignment.total_len>=min_match) {
					list_file<<read_id<<'\t'<<read_len<<'\t'<<alignment.read_start<<'\t'<<alignment.read_end<<'\t'
					        <<adpt_id<<'\t'<<adpt_len<<'\t'<<alignment.adpt_start<<'\t'<<alignment.adpt_end<<'\t'
					        <<alignment.total_len<<'\t'<<alignment.total_mis<<'\t'<<alignment.total_gap<<endl;
					adpt_num++;
					if (alignment.read_start<=3) {
						empty_num++;
					}
				}
			}

		}

		fastq_file.clear();
		fastq_file.seekg(0,ios::beg);

		float seq_percent = ((float) adpt_num)/read_num*100;
		float empty_percent = ((float) empty_num)/read_num*100;
		stat_file<<adpt_id<<'\t'<<adpt_num<<" ("<<setprecision(3)<<fixed<<seq_percent<<"%)"<<'\t'<<empty_num<<" ("<<empty_percent<<"%)"<<'\t'<<adpt<<endl;

	}

	stat_file<<"\ntotal_reads: "<<read_num<<endl;
	fastq_file.close();
	adpt_file.close();
	list_file.close();
	stat_file.close();
}



  /********************************************************/
 /*********		Sub Function		**********/
/********************************************************/
alignment get_alignment(string *read, int read_pos, string *adpt, int adpt_pos, int seed_len, int max_mis) {

	int adpt_len=(*adpt).length();
	int read_len=(*read).length();
	alignment alignment;
	alignment.total_gap=0;

	max_mis<0?0:max_mis;

	int l_mis_num=0,r_mis_num=0;
	int *l_mis,*r_mis;
	l_mis = new int [max_mis+2];
	r_mis = new int [max_mis+2];
	for (int i=0; i<max_mis+2; i++) {
		l_mis[i]=0;
		r_mis[i]=0;
	}
	for(int i=1;;i++) {//get the left mismatch position
		if (read_pos-i<0 || adpt_pos-i <0) {
			l_mis_num++;
			l_mis[l_mis_num]=i;
			break;
		}
		if ((*read)[read_pos-i] != (*adpt)[adpt_pos-i]) {
			l_mis_num++;
			l_mis[l_mis_num]=i;
		}
		if (l_mis_num>=max_mis+1) {
			break;
		}

	}

	for (;l_mis_num > 1 && l_mis[l_mis_num]==l_mis[l_mis_num-1]+1;) {//continuous mismatches at the end of the macth should be deleted;
		l_mis_num--;
	}

	for(int i=1;;i++) {//get the right mismatch position
		if (read_pos+seed_len+i > read_len || adpt_pos+seed_len+i > adpt_len) {
			r_mis_num++;
			r_mis[r_mis_num]=i;
			break;
		}
		if ((*read)[read_pos+seed_len-1+i] != (*adpt)[adpt_pos+seed_len-1+i]) {
			r_mis_num++;
			r_mis[r_mis_num]=i;
		}
		if (r_mis_num>=max_mis+1) {
			break;
		}
	}

	for (;r_mis_num > 1 && r_mis[r_mis_num]==r_mis[r_mis_num-1]+1;) {//continuous mismatches at the end of the macth should be deleted;
		r_mis_num--;
	}

/***********************************************************************************/

	int max_len=0,match_len=0;
	int l_mis_id,r_mis_id,l_max_id,r_max_id;

	if (l_mis_num+r_mis_num <= max_mis+2) {
		l_mis_id=l_mis_num;
		r_mis_id=r_mis_num;
		max_len = l_mis[l_mis_id]-1+r_mis[r_mis_id]-1;
		l_max_id=l_mis_id-1;
		r_max_id=r_mis_id-1;
	}else{
		for (int i=l_mis_num;i>=1;i--){
			l_mis_id = i;
			r_mis_id = max_mis+2-i;
			match_len=l_mis[l_mis_id]-1+r_mis[r_mis_id]-1;
			if (match_len>=max_len) {
			max_len=match_len;
			l_max_id=l_mis_id-1;
			r_max_id=r_mis_id-1;
			}
		}
	}

/***********************************************************************************/

	alignment.total_len=max_len+seed_len;
	alignment.total_mis = l_max_id + r_max_id;
	alignment.read_start = read_pos-(l_mis[l_max_id+1]-1);
	alignment.read_end = alignment.read_start + alignment.total_len - 1;
	alignment.adpt_start = adpt_pos-(l_mis[l_max_id+1]-1);
	alignment.adpt_end = alignment.adpt_start+alignment.total_len - 1;
	alignment.mis_rate =((float) alignment.total_mis) /((float) alignment.total_len);
	delete l_mis;
	delete r_mis;
	return (alignment);//match successfully

}

void usage() {
	cout<<"\nAuthor: BENM <binxiaofeng@gmail.com>\n"
		<<"Version: 0.1.2\n"
		<<"Date: 2012-11-22\n";
       cout << "\nUsage: scanAP -i <*.fq> -a <*.fa> -s <*.fq.stat> -d <*.fq.detail> [options]\n"
			<< "  -i <str>   input fastq file of reads\n"
			<< "  -a <str>   input fasta file of adapters\n"
			<< "  -s <str>   output statistics file\n"
			<< "  -d <str>   output alignment detail file\n"
			<< "  -k <int>   set length of seed, default:10\n"
			<< "  -e <int>   set the maximum mismatches, default:3\n"
			<< "  -m <int>   set the minimal match length, default:10\n"
			<< "  -l         open local alignment;\n"
			<< "  -h         output help information\n"
			<< endl ;
	cout<<"detail output(separated by tab):\n"
		<<"reads_id reads_len reads_start reads_end AP_id AP_len AP_start AP_end align_len mismatch gap"<<endl;
	cout<<"statfile(separated by tab):\n"
		<<"AP_id polluted_reads empty_reads AP_sequence\n"
		<<"coordinate start calculated from 0\n"<<endl;
	exit(1);
}

//get the score matrix and gap penalty from matrix file:./align.mat
void get_scoreMatrix(int &gapPenalty,map< char, map<char,int> > &scoreMatrix) {
	fstream mat_file("align.mat");
	string mat_line;
	vector<string> gapline,vecline0;
	gapPenalty = -5;//default gap penalty
// 	map< char, map<char,int> > scoreMatrix;
	if (mat_file==NULL)	{
		// cout<<"**Can not find the configure file:\"align.conf\"**"<<endl;
		// cout<<"**Check the directory of this program again**" <<endl;
		gapPenalty=-10;
		char base[]={'A','C','G','T','N'};
		int size=sizeof(base)/sizeof(*base);
		for (int i=0;i<size;i++){
			for (int j=0;j<size;j++){
				if (i==size-1 || j==size-1){
					scoreMatrix[base[i]][base[j]]=0;
				}else{
					if (i==j && j<size-1){
						scoreMatrix[base[i]][base[j]]=2;
					}else{
						if (abs(i-j)%2==0){
							scoreMatrix[base[i]][base[j]]=-5;
						}else{
							scoreMatrix[base[i]][base[j]]=-7;
						}
					}
				}
			}
		}
	}else{
		getline(mat_file,mat_line,'\n');
		split(mat_line,gapline,"\t");
		gapPenalty=atoi(gapline[1].c_str());//read the gap penalty;
		getline(mat_file,mat_line,'\n');//read the title of the scorematrix;
		split(mat_line,vecline0,"\t");
		while(getline(mat_file,mat_line,'\n'))	{//read sequences from two fastaq files;
			vector<string> vecline;
			split(mat_line,vecline,"\t");
			scoreMatrix[vecline[0][0]][vecline0[0][0]]=atoi(vecline[1].c_str());
			scoreMatrix[vecline[0][0]][vecline0[1][0]]=atoi(vecline[2].c_str());
			scoreMatrix[vecline[0][0]][vecline0[2][0]]=atoi(vecline[3].c_str());
			scoreMatrix[vecline[0][0]][vecline0[3][0]]=atoi(vecline[4].c_str());
			scoreMatrix[vecline[0][0]][vecline0[4][0]]=atoi(vecline[5].c_str());
		}
		mat_file.close();//close the configure file!
	}
}
/* ##gapPenalty	-10
*	A		C		G		T		N
*A	2		-7		-5		-7		0
*C	-7		2		-7		-5		0
*G	-5		-7		2		-7		0
*T	-7		-5		-7		2		0
*N	0		0		0		0		0
*/
