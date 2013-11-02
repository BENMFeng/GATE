#include<iostream>
#include<fstream>
#include<string>

using namespace std;

const char * README =
"argv1 input noclosure scaffold file\n"
"output to STDOUT\n";

void Output(string & name, string & seq)
{
	for(register size_t i(0),j(0); i < seq.size() && j < seq.size();)
	{
		if(seq[i] == 'N')
			i = seq.find_first_not_of("N", i);
		else
		{
			j = seq.find_first_of("N", i);
			cout<<name<<"_"<<i<<"_"<<((j>=seq.size())?seq.size():j)<<endl
				<<seq.substr(i, j - i)<<endl;
			i = j;
		}
	}
}

int main(int argc, char ** argv)
{
	if(argc != 2)
	{
		cerr<<README;
		return 1;
	}
	
	ifstream I(argv[1]);
	if(!I)
	{
		cerr<<"Error opening file: "<<argv[1]<<endl;
		return 2;
	}

	string scafname;
	string seq;
	string line;
	for(;;)
	{
		getline(I, line);
		if(!I)
		{
			Output(scafname,seq);
			break;
		}
		if(line.empty())
		{
			Output(scafname,seq);
			seq.clear();
			continue;
		}
		if(line[0]=='>')
		{
			Output(scafname, seq);
			seq.clear();
			scafname=line;
			continue;
		}
		seq += line;
	}
}
