/*file:checkfastaq.cpp
 //*test file:s_1_sequence.txt
 //*author: Huang Quanfei, email: huangqf@genomics.org.cn, date: 2008-3-21
        //  Fan Wei,       email: fanw@genomics.org.cn, date: 2008-3-21
        // Revised by: BENM email: binxiaofeng@gmail.com, date: 2013-3-22
 //*function: (1) caculate distribution of base vs cycle; (2) caculate distribution of quality vs cycle.
 // */


#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>
#define _A 0
#define _C 1
#define _G 2
#define _T 3
#define _N 4
#define MAX_LENGTH 1000

using namespace std;

void usage();

int main (int argc, char *argv[ ])
{
	fstream file1;
	int nseq=0;//number of sequences
	int length=0;//length of a sequence
	int qMAX=0;// max number of qualities
	int lengthmax=0;//max length of sequences
	int total=0;//total number of bases

	static int sum[5];

	static int qsum[256];
	static int pqsum[MAX_LENGTH][256];
	static int psum[MAX_LENGTH][5];
	double err=0.;
	double A_T=0.;
	double G_C=0.;
	int q10=0;
	int q20=0;
	int q30=0;
	int q40=0;
	double q_err[256];
	int ifrag = 0;
	string readfile;
	int qshift = 33;
	int opt;

	while ((opt=getopt(argc, argv, "i:q:h"))!=-1) {
		switch (opt) {
			case 'i': ifrag=1; readfile=optarg; break; //input file
			case 'q': qshift=atoi(optarg); break;
			case 'h': usage(); break;
			default:  usage();
		}
	}

	if (!(ifrag))	{
		cout<<"Please input the name of a fastaq file."<<endl;
		usage();
	}

	fstream fastq_file;
	fastq_file.open (readfile.c_str(), ios::in);

	if (fastq_file==NULL)	{
		cout<<"Fail to open fastaq file"<<endl;
		usage();
	}

	for(int i=0;i<256;i++){
		q_err[i]=1.0/(1.0 + pow((double)10.0,(double)((i-qshift)*0.1)));
	}


	string textline;
	while (getline(fastq_file,textline,'\n'))	{
		if (textline[0]=='@')	{
			int q_length=0;
			string seq, qual, id;
			getline(fastq_file,seq,'\n');
			getline(fastq_file,id,'\n');
			getline(fastq_file,qual,'\n');

/*		length of the sequence		*/
			length=seq.length();
			q_length=qual.length();

			if ((seq[0]=='A'|seq[0]=='C'|seq[0]=='G'|seq[0]=='T'| seq[0]=='N') && (length==q_length))
			{
				++nseq;//number of sequence
				if (length>lengthmax)
					lengthmax=length;

/*		sequence			*/
				for (int i=0;i<length;i++)	{
					switch (seq[i])	{
						case 'A':
							++psum[i][_A];break;
						case 'C':
							++psum[i][_C];break;
						case 'G':
							++psum[i][_G];break;
						case 'T':
							++psum[i][_T];break;
						case 'N':
							++psum[i][_N];break;
					}
				}
				total+=length;

/*			qualities			*/
				for (int i=0;i<length;i++)	{
					int q = qual[i]-qshift;
					if(q>0){
						if(q>qMAX)
							qMAX=q;
						++pqsum[i][q];
					}
					else{
						++pqsum[i][0];
					}

					if (seq[i]!='N')
					{
						if (q>=40) {++q40;++q20;++q10;}
						else if (q>30) {++q30;++q20;++q10;}
						else if (q>=20) {++q20;++q10;}
						else if (q>=10) {++q10;}
						err+=q_err[q+64];
					}
				}

			}

		}
	}
	fastq_file.close();
	for(int i=0;i<lengthmax;i++){
		sum[_A]+=psum[i][_A];
		sum[_C]+=psum[i][_C];
		sum[_G]+=psum[i][_G];
		sum[_T]+=psum[i][_T];
		sum[_N]+=psum[i][_N];
	}

	for(int i=0;i<lengthmax;i++) {
		for(int j=0;j<=qMAX;j++) {
			qsum[j]+=pqsum[i][j];
		}
	}

	 A_T = (100.0*(sum[_A]-sum[_T]))/total;
	 G_C = (100.0*(sum[_G]-sum[_C]))/total;

/*			output			*/
 	if(total)	{
		cout<<" the default quality shift value is: "<<qshift<<", "
			<<nseq<<" sequences, "
			<<total<<" total length, Max length:"
			<<lengthmax<<", average length:"
			<<setiosflags(ios::fixed)<<setprecision(2)<<(double)total/nseq<<endl;
		cout<<"Standard deviations at 0.25:  total "<<100*(sqrt(0.25*(double)total)/total)
			<<"%, per base "<<100*(sqrt(0.25*(double)nseq)/nseq)
			<<"%"<<endl;
		cout<<"             A     C     G     T     N ";
		for(int i=0;i<=qMAX;i++)
			cout<<setw(4)<<i<<' ';
		cout<<endl;
		cout<<"Total    ";
		for(int i=_A;i<=_N;i++)
			cout<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<100*(double)sum[i]/total<<' ';
		for(int i=0 ;i<=qMAX;i++)
			cout<<setw(4)<<(int)(1000*((double)qsum[i]/total))<<" ";
		cout<<endl;
		for(int i=0;i<=lengthmax-1;i++)	{
			cout<<"base "<<setw(3)<<i+1<<' ';
			for(int j=_A;j<=_N;j++)
				cout<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<100*(double)psum[i][j]/nseq<<' ';
			for(int k=0;k<=qMAX;k++)
				cout<<setw(4)<<(int)(1000*((double)pqsum[i][k]/nseq))<<' ';
			cout<<endl;
		}
		cout<<endl;
		cout<<"%GC\tQ10\tQ20\tQ30\tQ40\t%A-%T\t%G-%C\tError Rate"<<endl;
		cout.precision(2);
		cout<<(100.0*(sum[_C]+sum[_G]))/(total-sum[_N])<<"\t"<<(100.0*q10)/total<<"\t"<<(100.0*q20)/total<<"\t"<<(100.0*q30)/total<<"\t"<<(100.0*q40)/total<<"\t"<<A_T<<"\t"<<G_C<<"\t"<<(100.0*err)/total<<endl;
	}

	return(0);
}

void usage() {
	cout<<"\nAuthor: BENM <binxiaofeng@gmail.com>\n"
		<<"Version: 0.1.0\n"
		<<"Date: 2013-03-22\n";
       cout << "\nUsage: checkfastq -i <*.fq> [-f L]\n"
			<< "  -i <file>	input fastq file of reads\n"
			<< "  -q <int>	FASTQ quality shift, default: 33;\n"
			<< "  -h		output help information\n"
			<< endl ;
	exit(1);
}
