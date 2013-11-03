#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <math.h>
using namespace std;
#define PI   3.1415926535897932384626433832795
#define E 2.7182818284590452353602874713527


const double a2[]={ 1.0/12.0, -1.0/360, 1.0/1260.0 };
double Fact(int n);
double P_function(int N1,int N2,int x,int y);
double P_value_cacul(int N1,int N2,int x,int y);
void usage();
int opt;

int main(int argc, char* argv[])
{
      if(argc == 1)
		usage();
	else{
        	char  *inputfile = 0;
        	char  *outfile = 0;
        	int N1;
        	int N2;
		  while((opt = getopt(argc, argv, "i:o:")) != -1 )
	  		{

		 		switch(opt)
		  		{
			  	case 'i':
					inputfile = optarg;
					break;
            		  	case 'o':
               				outfile = optarg;
					break;
			  	case '?':
					printf("unknown option: %c\n", optopt);
				 	return EXIT_SUCCESS;
					break;
		  		}
	  	 }//end while opt

          ifstream INf (inputfile);
          ofstream Fout(outfile);

       if(!INf)
        {
            cerr <<"unable to open file"<<INf<<"\n";
            exit (-1);
        }
        else
        {
            string buf;

		
            while(!INf.eof())
            {
               getline(INf,buf);
               if(buf!="")
		{
		   string  col[5];	
		   int pre_pos=0;
		   int next_pos=buf.find_first_of("\t",0);
		   for(int i=0;i<4;i++)
		     {
			 col[i]=buf.substr(pre_pos,next_pos);
			 pre_pos=next_pos;
			 next_pos=buf.find_first_of("\t",pre_pos+1);

			}
			col[4]=buf.substr(pre_pos,buf.size()-1);
		   
			//col[4]=buf.substr(pre_pos,next_pos);	
               	
		   int num[4];
	    	   for(int i=0;i<4;i++)
		   {
			num[i]=atoi(col[i+1].c_str());
			//cout<<num[i]<<endl;
		   }
		   double p_value=P_value_cacul(num[0],num[1],num[2],num[3]);
		  string up_down;
		 if(p_value>0.5)
		 {p_value=1-p_value;
		  up_down="+";	}
		  else
		 {
 			up_down="-";
		 }
		 if(p_value<0){p_value=0;}
               	   Fout<<col[0]<<"\t"<<num[0]<<"\t"<<num[1]<<"\t"<<num[2]<<"\t"<<num[3]<<"\t"<<p_value<<"\t"<<up_down<<endl;
		  	 
		}
            }


       }

     Fout.close();
     INf.close();

    }

 return 0;
}

double P_value_cacul(int N1,int N2,int x,int y)
{
 double p_value=0;
 
 for (int i=0;i<=y;i++)
    {
        double temp_p=P_function(N1,N2,x,i);
        p_value+=pow(10,temp_p);
    }
 

 return p_value;
}


// log10((N2/N1)^y*((x+y)!/x!*y!)(1+N2/N1)^(x+y+1))
//= y*log(N2/N1) +log((x+y)!)-log(x!)-log(y!)-(x+y+1)log[1+(N2/N1)] ;
double P_function(int N1,int N2,int x,int y)
{

  double part1,part2,part3,part4,part5;
  part1=y*(log10(double(N2)/double(N1)));
  part2=Fact(x+y);
  part3=Fact(x);
  part4=Fact(y);
  part5=(x+y+1)*log10(1+double(N2)/double(N1));
  double result=part1+part2-part3-part4-part5;
 
return result;
}



double Fact(int n)
{
     double factorial;
    if(n!=0)
    {
       double logR;
       double s, //sum of a2
       item;       //each item of a2
       int i;
	logR=0.5*log(2.0*PI)+((double)n+0.5)*log(n)-(double)n;

       //¡¡s= (1/12/n -1/360/n^3 + 1/1260/n^5)
       for (item=1/(double)n,s=0.0,i=0;i<sizeof(a2)/sizeof(double);i++)
       {
              s+= item * a2[i];
              item /= (double)(n)* (double)n; //item= 1/(n^(2i+1))
       }
       logR+=s;

       //log10(n)=ln(n)/ln(10)
      factorial=logR/log(10);

    }
    else
    {factorial=0;
	//0!=1,log10(1)=0;
      }

    return factorial;
}

void usage()
{

  cout << "Calculate the significance of each gene expression in two sample " <<endl;
  cout << "Options:" << endl;
  cout << "\t\t-i STR: inputfile " <<endl;
  cout << "\t\t\tinput formate:  genen_id \\t smaple1_num \\t sample2_num \\t express_num1 \\t express_num2 " <<endl;
  cout << "\t\t-o STR: outputfile" <<endl;
   
}

