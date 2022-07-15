// FOGSAA algorithm implementation with affine gap penalty
//               Developed  By Angana Chakraborty



// -------------------------------------------
 // ############      DNA/Gene Sequences    ###############
/*-----------------------------------------
Command line parameters (for alignments without affine gap penalty)
./fogsaa f1 f2 d 0 mt mst gp
---- where f1 --> name of the 1st file.
           f2 --> name of the 2nd file
           d  --> 1 for DNA sequences 
           0  --> as without affine gap penalty  
           mt --> score of a match
           mst--> score for mismatch
           gp --> gap penalty
    example ./fogsaa seq1.txt seq2.txt 1 0 1 -1 -2
-----------------------------------------------------------
-----------------------------------------------------------
-----------------------------------------------------------
Command line parameters ( for alignments with affine gap penalty)
./fogsaa f1 f2 d 1 mt mst go ge 
---- where f1 --> name of the 1st file.
           f2 --> name of the 2nd file
           d  --> 1 for DNA sequences 
           1  --> as with affine gap penalty 
           mt --> score of a match
           mst--> score for mismatch
           go --> gap open penalty
           ge --> gap extension penalty, Total penalty for a gap of length L= (go+L*ge)
    example ./fogsaa seq1.txt seq2.txt 1 1 1 -1 -10 -1
--------------------------------------------------------------------
--------------------------------------------------------------------

###################    Protein Sequences    #######################
-------------------------------------------------------------------
Command line parameters (for alignments without affine gap penalty)
./fogsaa f1 f2 d 0 gp
---- where f1 --> name of the 1st file.
           f2 --> name of the 2nd file
           d  --> 2 for Protein sequences 
           0  --> as without affine gap penalty 
           gp --> gap penalty
    example ./fogsaa seq1.txt seq2.txt 2 0 -2
    protein score matrix is BLOSUM64 stored on file fscore.txt
-------------------------------------------------------------------
-------------------------------------------------------------------
-------------------------------------------------------------------
Command line parameters ( for alignments with affine gap penalty)
./fogsaa f1 f2 d 1 go ge 
---- where f1 --> name of the 1st file.
           f2 --> name of the 2nd file
           d  --> 2 for Protein sequences 
           1  --> as with affine gap penalty 
           go --> gap open penalty
           ge --> gap extension penalty, Total penalty for a gap of length L= (go+L*ge)
    example ./fogsaa seq1.txt seq2.txt 2 1 -10 -1
    protein score matrix is BLOSUM64 (default) stored on file fscore.txt
---------------------------------------------------------------------------
###########################################################################
###########################################################################
*/


#include<iostream>
#include<stdio.h>
#include<string.h> 
#include<stdlib.h>

using namespace std;
#include <sys/time.h>
#include <unistd.h>






int m,n, expanded=0,score=0,lowerBound,a[3],b[3],ch[3],optp1,optp2,approximate=0,threshold, DNA, affine=0, Mt,Mst,Gp,Go,Ge; // a=lower b=upper c=child index
 // Mt == match score, Mst== Mismatch score, Gp== Gap penalty,, Go== Gap open penalty,, Ge== Gap extension penalty, 

// structure definition
int fscore[20][20],prmax,prmin;
char amino[25];
typedef struct celltype
	{
		int present_score,lower,upper,type,filled; // type '1' is for (i+1,j+1),, type '2' is for (i,j+1), type '4' is for (i+1,j)  to identify them distinctly
	}cell;

cell **c;


typedef struct qtype
{
	int p1,p2,type_upto_next,next_type,next_lower;
	struct qtype *next;
}queue;
int maxPointer=-1;
queue **q=NULL;

FILE *fp1,*fp2,*fp3,*fp4;
char * pair1,*pair2;
int curp1=0,curp2=0;

// function prototype
int ins_queue(int,int,int ,int,int ,int);
void del_queue(void);
void align(int,int);
void sort(void);
void calculate_score(int p,int *l,int *u,int p1,int p2);

int main(int argc, char * argv[])
{	char ca,cb,pathend=0,loop=1;
 int type_total=1,new_type,new_score,np1,np2,new_lower,new_upper,next_lower,next_upper;
 int p,present;
 	
	
	// parameter checking
	if((strcmp(argv[3],"1")==0) && (strcmp(argv[4],"0")==0))  // gene seq without affine gap
	{  
	   if(argc!=8)
	    printf(" Invalid argument\n");
	    else
	    {
		DNA=1;
		affine=0;
	   	 Mt=atoi(argv[5]);
	   	 Mst=atoi(argv[6]);
	   	 Gp=atoi(argv[7]);
	    }
	}
	 else if((strcmp(argv[3],"1")==0) && (strcmp(argv[4],"1")==0)) // gene seq with affine gap
	 {
	 	if(argc!=9)
	 	printf(" Invalid argument\n");
	 	else
	 	{
	 		DNA=1;
			affine=1;
			Mt=atoi(argv[5]);
	 		Mst=atoi(argv[6]);
	 		Go=atoi(argv[7]);
	 		Ge=atoi(argv[8]);
	 	}
	 }
	 else if((strcmp(argv[3],"2")==0) && (strcmp(argv[4],"0")==0))  // protein seq without affine gap
	 {
	 	
 	
	 	if(argc!=6)
	 	printf(" Invalid argument\n");
	 	else
	 	{
	 		DNA=2;
			affine=0;
			Gp=atoi(argv[5]);
	 		
	 	}
	 }
	  else if((strcmp(argv[3],"2")==0) && (strcmp(argv[4],"1")==0))  // protein seq with affine gap
	  {	
	  	if(argc!=7)
	  	printf(" Invalid argument\n");
	  	else
	  	{	
	  		DNA=2;
			affine=1;
			Go=atoi(argv[5]);
	  		Ge=atoi(argv[6]);
	  	}
	  }
	  else
	  printf(" Invalid argument\n");
	  
 // loop for multiple pairs of sequences
 for(loop=1;loop<=1;loop++)
 { 
 

	fp1=fopen(argv[1],"r");
	fp2=fopen("seq1.txt","w");
	while((ca=getc(fp1))!='\n')
	{ }

	while((ca=getc(fp1))!=EOF)
	{
		if((ca>='A') &&(ca<='Z'))
			putc(ca,fp2);
		else if((ca>='a') &&(ca<='z'))
			putc(ca-32,fp2);
	}
	fclose(fp1);
	fclose(fp2);

	fp1=fopen(argv[2],"r");
	fp2=fopen("seq2.txt","w");
	while((ca=getc(fp1))!='\n')
	{ }

	while((ca=getc(fp1))!=EOF)
	{
		if((ca>='A') &&(ca<='Z'))
			putc(ca,fp2);
		else if((ca>='a') &&(ca<='z'))
			putc(ca-32,fp2);
	}
	fclose(fp1);
	fclose(fp2);
	// create the array of pairs
	
	
	int i;
	fp1=fopen("seq1.txt","r");
	fseek(fp1, 0L, SEEK_END);
	m = ftell(fp1);
	pair1 = (char *)malloc((m+1) * sizeof(char));
	fseek(fp1, 0L, SEEK_SET);
	cout<<m<<"\n";
	for(i=0;i<m;i++)
	{
		fseek(fp1,(long) (i), SEEK_SET);
		pair1[i]=getc(fp1);
		
		
	}
	pair1[i]='\0';
	fp2=fopen("seq2.txt","r");
	fseek(fp2, 0L, SEEK_END);
	n = ftell(fp2);
	pair2 = (char *)malloc((n+1) * sizeof(char));
	fseek(fp2, 0L, SEEK_SET);
	cout<<n<<"\n";
	for(i=0;i<n;i++)
	{
		fseek(fp2,(long) (i), SEEK_SET);
		pair2[i]=getc(fp2);
	}
	pair2[i]='\0';
	fclose(fp1);
	fclose(fp2);
	
	//for protein
	if(DNA==2)
	{

		// initialize the amino acid matrix
		fp1=fopen("amino.txt","r");
		char cc;

		for(int i=0;i<20;i++)
		{
			fscanf(fp1,"%c",&cc);
			
			amino[cc-'A']=i;
			
		}

		
		fclose(fp1);
		
		// read the score matrix from file
		fp1=fopen("fscore.txt","r");
		fscanf(fp1,"%d",&fscore[0][0]);
		prmin=fscore[0][0];
		prmax=fscore[0][0];
              	fseek(fp1, 0L, SEEK_SET);
                for(int i=0;i<20;i++)
		{
			
			for(int j=0;j<=i;j++)
			{
			 fscanf(fp1,"%d",&fscore[i][j]);
			 fscore[j][i]=fscore[i][j];
			 if(fscore[i][j]>prmax)
			 	prmax=fscore[i][j];
			 if(fscore[i][j]<prmin)
			 	prmin=fscore[i][j];
			 //printf("%d  ",fscore[i][j]);
			} 
			//printf("\n");
			 
		}
		//printf("max =%d and min =%d\n",prmax,prmin);
		
		fclose(fp1);
		
	}
		
	
	
	
	
	//get time
	struct timeval start, end;

    long mtime, seconds, useconds;    

    gettimeofday(&start, NULL);
    
    

	

	 c = (cell **)malloc((m+1) * sizeof(cell *));
	 for (i = 0; i <=m; i++) {
    c[i] =(cell *) malloc( (n+1)* sizeof(cell));
  }
		
	
	
	
	
	c[0][0].present_score=0;
	c[0][0].type=-1;
	
	calculate_score(c[0][0].present_score,&(c[0][0].lower),&(c[0][0].upper),0,0);
	lowerBound=c[0][0].lower;
	

	// estimate the threshold value with 30% similarity requirement
	int ml,sl,th;
	if(m>n)
	{
		ml=m;
		sl=n;
	}
	else
	{
		ml=n;
		sl=m;
	}
	th=ml*30/100;
	if(DNA==1)
	{
	if(affine==0)
	threshold=th*Mt+(sl-th)*Mst+Gp*(ml-sl);
	else
	threshold=th*Mt+(sl-th)*Mst+Ge*(ml-sl)+Go;
	}
	else  // protein
	{
	if(affine==0)
	threshold=th*prmax+(sl-th)*prmin+Gp*(ml-sl);
	else
	threshold=th*prmax+(sl-th)*prmin+Ge*(ml-sl)+Go;
	}
		
	
	

	int v=c[0][0].upper-c[0][0].lower+1;
	q=new queue*[v];  // creating array of pointers to queuetype
	for(i=0;i<v;i++)
		q[i]=NULL;




	if(m!=0 && n!=0)
	{
		do
		{
			
			pathend=1;
		while((curp1<=m-1) || (curp2 <=n-1))
		{
			
			present=c[curp1][curp2].present_score;
			// expand the child node
			
			if((type_total==1)||(type_total==2)||(type_total==4))
			
			{ 
				//this is the first child
				
			if(curp1<=m-1 && curp2<=n-1)
			{
					
					ca=pair1[curp1];
					cb=pair2[curp2];
					if(DNA==2)
					{
						
						Mt=fscore[amino[ca-'A']][amino[cb-'A']];
						//printf("The value of %c,%c at %d,%d=%d\n",ca,cb,amino[ca-'A'],amino[cb-'A'],Mt);
						Mst=Mt;	
					}
					if(ca==cb)
						p=Mt;
					else
						p=Mst;

			                 
					
					calculate_score(present+p,&a[0],&b[0],curp1+1,curp2+1);
					if(affine==0)
					{
					
						calculate_score(present+Gp,&a[1],&b[1],curp1,curp2+1);
						calculate_score(present+Gp,&a[2],&b[2],curp1+1,curp2);
					}
					else // affine gap penalty
					{
						if((c[curp1][curp2].type==1)||(c[curp1][curp2].type==-1))
						{
						calculate_score(present+Go+Ge,&a[1],&b[1],curp1,curp2+1);
						calculate_score(present+Go+Ge,&a[2],&b[2],curp1+1,curp2);
						}
						else if(c[curp1][curp2].type==2)// as the gap is already opened in the first chain
						{
						calculate_score(present+Ge,&a[1],&b[1],curp1,curp2+1);
						calculate_score(present+Go+Ge,&a[2],&b[2],curp1+1,curp2);
						}
						else // as the gap is already opened in the 2nd chain
						{
						calculate_score(present+Go+Ge,&a[1],&b[1],curp1,curp2+1);
						calculate_score(present+Ge,&a[2],&b[2],curp1+1,curp2);
						}
						
					}
					
					
			        sort(); // sort the children according to upper and when there is tie... check the lower
					
					
					
					if(ch[0]==1)
					{
						
						// child would be of type 1
						new_type=1;
						np1=curp1+1;
						np2=curp2+1;
						new_score=present+p;
						
						

					}
					else if(ch[0]==2)
					{	
						
						// child would be of type 2
						new_type=2;
						np1=curp1;
						np2=curp2+1;
						if(affine==0)
						new_score=present+Gp;
						else // in the affine gap penalty scheme
						{
							if(c[curp1][curp2].type==1 || c[curp1][curp2].type==4 ||(c[curp1][curp2].type==-1))  // the child is opening a new gap now
						
								new_score=present+Go+Ge;
							else if(c[curp1][curp2].type==2)// the gap was already opened in first chain
								new_score=present+Ge;
							
							
						}
						
						


						
					}
					else
						{
						// child would be of type 4
						
						new_type=4;
						np1=curp1+1;
						np2=curp2;
						if(affine==0)
						new_score=present+Gp;
						else // in the affine gap penalty scheme
						{
							if(c[curp1][curp2].type==1 || c[curp1][curp2].type==2 || (c[curp1][curp2].type==-1))  // the child is opening a new gap now in 1st chain
						
								new_score=present+Go+Ge;
							else // the gap was already opened in 2nd chain
								new_score=present+Ge;
							
							
						}
						
						

						
					}

					//printf("before queue\n");
					maxPointer=ins_queue(curp1,curp2,new_type+ch[1],ch[1],a[1],b[1]);
					//printf("after queue\n");
			}
			else if(curp1<=m-1)
			{
				// only type '4' child is possible
						new_type=4;
						np1=curp1+1;
						np2=curp2;
						if(affine==0)
						new_score=present+Gp;
						else // in the affine gap penalty scheme
						{
							if(c[curp1][curp2].type==1 || c[curp1][curp2].type==2 || (c[curp1][curp2].type==-1))  // the child is opening a new gap now
						
								new_score=present+Go+Ge;
							else // the gap was already opened
								new_score=present+Ge;
							
							
						}
						
				
			}
			else
			{
				// only type '2' child is possible
						new_type=2;
						np1=curp1;
						np2=curp2+1;
						if(affine==0)
						new_score=present+Gp;
						else // in the affine gap penalty scheme
						{
							if(c[curp1][curp2].type==1 || c[curp1][curp2].type==4 || (c[curp1][curp2].type==-1))  // the child is opening a new gap now
						
								new_score=present+Go+Ge;
							else // the gap was already opened
								new_score=present+Ge;
							
							
						}
						
				

			}
			
			
			}  // end of first child
			else if((type_total==3)||(type_total==5)||(type_total==6))
			{
				// this is the 2nd child
				
				if(new_type==1)
				{
					// 2nd child is of type type 1
					np1=curp1+1;
					np2=curp2+1;
					ca=pair1[curp1];
					cb=pair2[curp2];
					if(DNA==2)
					{
						Mt=fscore[amino[ca-'A']][amino[cb-'A']];
						//printf("The value =%d\n",Mt);
						Mst=Mt;	
					}
					if(ca==cb)
						p=Mt;
					else
						p=Mst;
					new_score=present+p;
					if(7-type_total==2)
					{
						if(affine==0)
						calculate_score(present+Gp,&next_lower,&next_upper,curp1,curp2+1);
						else
						{
							if(c[curp1][curp2].type==1 || c[curp1][curp2].type==4 || (c[curp1][curp2].type==-1))
							 calculate_score(present+Go+Ge,&next_lower,&next_upper,curp1,curp2+1);
							 else
							 calculate_score(present+Ge,&next_lower,&next_upper,curp1,curp2+1);
						}
					}
					else if(7-type_total==4)
					{
						if(affine==0)
						calculate_score(present+Gp,&next_lower,&next_upper,curp1+1,curp2);
						else
						
						{
							if(c[curp1][curp2].type==1 || c[curp1][curp2].type==2 || (c[curp1][curp2].type==-1))
							calculate_score(present+Go+Ge,&next_lower,&next_upper,curp1+1,curp2);
							else
							calculate_score(present+Ge,&next_lower,&next_upper,curp1+1,curp2);
							
						}
					}
					//printf("before queue2\n");
					maxPointer=ins_queue(curp1,curp2,7,7-type_total,next_lower,next_upper);
					//printf("before queue2\n");

				}
				else if(new_type==2)
				{
					//2nd child is of type type 2
					np1=curp1;
					np2=curp2+1;
					if(affine==0)
					new_score=present+Gp;
					else
					{	
						if(c[curp1][curp2].type==1 || c[curp1][curp2].type==4 || (c[curp1][curp2].type==-1))
						   new_score=present+Go+Ge;
						else 
						   new_score=present+Ge;  // because gap has already been opened
							
					}
					if(7-type_total==1)
					{
						ca=pair1[curp1];
						cb=pair2[curp2];

						if(DNA==2)
						{
							Mt=fscore[amino[ca-'A']][amino[cb-'A']];
							//printf("The value =%d\n",Mt);
							Mst=Mt;	
						}
						if(ca==cb)
							p=Mt;
						else
							p=Mst;
						calculate_score(present+p,&next_lower,&next_upper,curp1+1,curp2+1);
					}
					else if(7-type_total==4)
					{
					  if(affine==0)
					    calculate_score(present+Gp,&next_lower,&next_upper,curp1+1,curp2);
					  else
					    {
						    if(c[curp1][curp2].type==1 || c[curp1][curp2].type==2 || (c[curp1][curp2].type==-1))
						calculate_score(present+Go+Ge,&next_lower,&next_upper,curp1+1,curp2);
						else
						calculate_score(present+Ge,&next_lower,&next_upper,curp1+1,curp2);
					    }
					}
					//printf("before queue3\n");
					maxPointer=ins_queue(curp1,curp2,7,7-type_total,next_lower,next_upper);
					//printf("before queue3\n");
					
				}
				else
				{
					// 2nd child is of type type 4
					np1=curp1+1;
					np2=curp2;
					if(affine==0)
					new_score=present+Gp;
					else
					{	
						if(c[curp1][curp2].type==1 || c[curp1][curp2].type==2 || (c[curp1][curp2].type==-1))
						   new_score=present+Go+Ge;
						else 
						   new_score=present+Ge;  // because gap has already been opened
							
					}
					if(7-type_total==1)
					{
						ca=pair1[curp1];
						cb=pair2[curp2];
						if(DNA==2)
						{
							Mt=fscore[amino[ca-'A']][amino[cb-'A']];
							Mst=Mt;	
						}
						if(ca==cb)
							p=Mt;
						else
							p=Mst;
						calculate_score(present+p,&next_lower,&next_upper,curp1+1,curp2+1);
					}
					else if(7-type_total==2)
					{
					  if(affine==0)
					    calculate_score(present+Gp,&next_lower,&next_upper,curp1,curp2+1);
					  else
					    {
						    if(c[curp1][curp2].type==1 || c[curp1][curp2].type==4 || (c[curp1][curp2].type==-1))
						calculate_score(present+Go+Ge,&next_lower,&next_upper,curp1,curp2+1);
						else
						calculate_score(present+Ge,&next_lower,&next_upper,curp1,curp2+1);
					    }
					}
					//printf("before queue4\n");
					maxPointer=ins_queue(curp1,curp2,7,7-type_total,next_lower,next_upper);
					//printf("before queue4\n");

				}



			}// end of 2nd child
			else if(type_total==7)
			{
				// this is the third child
				

				if(new_type==1)
				{
					np1=curp1+1;
					np2=curp2+1;
					ca=pair1[curp1];
					cb=pair2[curp2];
					if(DNA==2)
					{
						Mt=fscore[amino[ca-'A']][amino[cb-'A']];
						//printf("The value =%d\n",Mt);
						Mst=Mt;	
					}
					if(ca==cb)
						p=Mt;
					else
						p=Mst;
					new_score=present+p;
				}
				else if(new_type==2)
				{
					np1=curp1;
					np2=curp2+1;
					if(affine==0)
					new_score=present+Gp;
					else
					{	
						if(c[curp1][curp2].type==1 || c[curp1][curp2].type==4 || (c[curp1][curp2].type==-1))
						   new_score=present+Go+Ge;
						else 
						   new_score=present+Ge;  // because gap has already been opened
							
					}
				}
				else
				{
					np1=curp1+1;
					np2=curp2;
					if(affine==0)
					new_score=present+Gp;
					else
					{	
						if(c[curp1][curp2].type==1 || c[curp1][curp2].type==2 || (c[curp1][curp2].type==-1))
						   new_score=present+Go+Ge;
						else 
						   new_score=present+Ge;  // because gap has already been opened
							
					}
				}

				// here no node is to be inserted in the queue
			}
			

			
			// write the new child node
			
			if(c[np1][np2].type<=4 && c[np1][np2].present_score>=new_score && c[np1][np2].filled==1) // skip the path if already expanded (same)node is better
			{
				 // printf("pruned already better at(%d,%d) old=%d,new=%d\n",np1,np2,c[np1][np2].present_score,new_score);
					pathend=0;
					
					break;
			}
		
			else  // insert in the cell
			{
				calculate_score(new_score,&new_lower,&new_upper,np1,np2);
				c[np1][np2].present_score=new_score;
				c[np1][np2].lower=new_lower;
				c[np1][np2].upper=new_upper;
				c[np1][np2].type=new_type;
				c[np1][np2].filled=1;
				
			}

		
			// print the node
				
		
			/*char c1,c2;
			
				ca=pair1[np1-1];
					cb=pair2[np2-1];
			if(np1==0)
				c1=' ';
			else if(curp1==np1)
				c1='-';
			else
				c1=ca;

			if(np2==0)
				c2=' ';
			else if(curp2==np2)
				c2='-';
			else
				c2=cb;
			printf("node(%c,%c,%d,%d,%d,%d,%d)\n",c1,c2,np1,np2,c[np1][np2].present_score,c[np1][np2].lower,c[np1][np2].upper);*/
			
			
			// now point the child as the current node
				curp1=np1;
				curp2=np2;
				type_total=1;

			if(c[np1][np2].upper<lowerBound)  
			{
				pathend=0;
				break;
			}
			
			
			expanded++; // counts no of nodes expanded


			
      

		}// end of while

		

	if((c[curp1][curp2].present_score>lowerBound) && pathend==1) // if the current path is not totally expanded , don't reset the lowerbound
	{
		lowerBound=c[curp1][curp2].present_score;
		optp1=curp1;
		optp2=curp2;
		//printf("Path end : new lowerBound=%d\n",lowerBound);
	}
	
	/*if((expanded>m*n/3) && (lowerBound<threshold))
			{
				printf("The given sequences are not globally similar(below 30%% of similarity)-- try local allignment --\n");
				printf(" The score below is near optimal\n");
				printf("The threshold== %d\n",threshold);
				approximate=1;
				printf("now the future is= %d\n",new_upper);
				break;
			}*/


	
	
			if(maxPointer!=-1 )
			{
				curp1=q[maxPointer]->p1;
				curp2=q[maxPointer]->p2;
				
				type_total=q[maxPointer]->type_upto_next;
				new_lower=q[maxPointer]->next_lower;
				new_upper=maxPointer+c[0][0].lower;
				new_type=q[maxPointer]->next_type;
				del_queue();
				
				
			}
			//printf(" The future is=%d\n",new_upper);
			int cpl;
			if(curp1>curp2)
			
				cpl=curp1;
				
			else
				cpl=curp2;
				
		if(((cpl>70*ml/100)&&(lowerBound<threshold))||(new_upper<threshold))
			{
				printf("The given sequences are not globally similar(below 30%% of similarity)-- try local allignment --\n");
				printf(" The score below is near optimal\n");
				printf("The threshold== %d\n",threshold);
				approximate=1;
				printf("now the future is= %d\n",new_upper);
				break;
			}
		
		}while(lowerBound<new_upper);


		// end timer
		gettimeofday(&end, NULL);

    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;

    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

    printf("Elapsed time: %ld milliseconds\n", mtime);


		// print the string
		fp3=fopen("new_alseq1.txt","w");
		fp4=fopen("new_alseq2.txt","w");
		
			align(optp1,optp2);
			
			fclose(fp3);
			fclose(fp4);

			printf("total nodes expanded==%d\n\n",expanded);
			
			printf("score= %d\n",c[optp1][optp2].present_score);
	}// end of if
	
	
	} // end of multiple pair loop
	
}// end of main



void align(int p1,int p2)
{
	
	char ca,cb,c1,c2;
	int pp1,pp2;
		ca=pair1[p1-1];
		cb=pair2[p2-1];
		
		
	if(c[p1][p2].type==-1)
	{
		
		return;
	}
	else
	{
		if(c[p1][p2].type==1)
		{
			pp1=p1-1;
			pp2=p2-1;
		}
			
		else if(c[p1][p2].type==2)
		{
			pp1=p1;
			pp2=p2-1;
		}
		else
		{
			pp1=p1-1;
			pp2=p2;
		}

		align(pp1,pp2);
		if(p1==m+1)
		{
			putc('-',fp3);
			c1='-';
		}
		else
		{
			if(pp1==p1)
			{
				putc('-',fp3);
				c1='-';
			}
			else
			{
				putc(ca,fp3);
				c1=ca;
			}
		}

		if(p2==n+1)
		{
			putc('-',fp4);
			c2='-';
		}
		else
		{
			if(pp2==p2)
			{
				putc('-',fp4);
				c2='-';
			}
			else
			{
				putc(cb,fp4);
				c2=cb;
			}
		}
		
      
	} // end of else

}


void del_queue()
{
	queue * temp;
	if(maxPointer==-1)
		return;
	else
	{
		temp=q[maxPointer];
		q[maxPointer]=temp->next;
		free(temp);
		while(q[maxPointer]==NULL)  // this row is empty , point next row
			maxPointer--;
	}
}


int ins_queue(int p1,int p2,int type_total,int next_type,int next_lower,int next_upper)
{
	int inserted=0;
	
	queue *p,*prev,*newnode;
	if((next_upper-c[0][0].lower)>=0)
	{
	newnode=(queue*)malloc(sizeof(queue));
	newnode->p1=p1;
	newnode->next_type=next_type;
	newnode->next_lower=next_lower;
	newnode->type_upto_next=type_total;
	newnode->p2=p2;
	//printf("queue insertion of =(%d,%d) type %d at %d\n",p1,p2,next_type,next_upper);
	
	if(maxPointer==-1)
	{
		q[next_upper-c[0][0].lower]=newnode;
		newnode->next=NULL;
		maxPointer=next_upper-c[0][0].lower;
	}
	else
	{
		if(q[next_upper-c[0][0].lower]==NULL)  // first memeber in this row
		{
			q[next_upper-c[0][0].lower]=newnode;
			newnode->next=NULL;
			if((next_upper-c[0][0].lower)>maxPointer)
					maxPointer=next_upper-c[0][0].lower;
		}
		else
		{
			// search the appropriate position in the row comparing the lower value
			p=q[next_upper-c[0][0].lower];
			
			prev=NULL;
			while(p!=NULL)
			{
				if(p->next_lower<=next_lower)
				{
					// insert here(before p)
					if(prev==NULL) // insert as the first row in the list
					{
						q[next_upper-c[0][0].lower]=newnode;
						newnode->next=p;
					}
					else
					{
						// insert in the middle of the row
						prev->next=newnode;
						newnode->next=p;
					}
					inserted=1;
					break;


				}
				else
				{
					prev=p;
					p=p->next;
				}
			}

			if(inserted==0)
			{
				// insert at the end
				prev->next=newnode;
				newnode->next=NULL;
			}
		}
	

	}
	}
	//else printf("skipped************************\n");
	return(maxPointer);
	
} // end of the function

void calculate_score(int p,int* l,int * u,int p1,int p2)
{
	if(DNA==2)
	{
		Mt=prmax;
		Mst=prmin;
	}
	if(affine==1)
	{
		Gp=Go+Ge;
	}
	if((m-p1)<=(n-p2))
	{
		*l=p+((m-p1)*Mst+Gp*((n-p2)-(m-p1)));
		*u=p+((m-p1)*Mt+Gp*((n-p2)-(m-p1)));
	}
	else
	{
		*l=p+Mst*(n-p2)+Gp*((m-p1)-(n-p2));
		*u=p+(n-p2)*Mt+Gp*((m-p1)-(n-p2));
	}
}



	

void sort(void)  // descending order
{
	int i,j,t;
	ch[0]=1;
	ch[1]=2;
	ch[2]=4;
	for(i=0;i<2;i++)
	{
		for(j=0;j<2-i;j++)
		{
			if((a[j]<a[j+1])||((a[j]==a[j+1])&&(b[j]<b[j+1])))
			{
				//swap
				t=a[j];
				a[j]=a[j+1];
				a[j+1]=t;

				t=ch[j];
				ch[j]=ch[j+1];
				ch[j+1]=t;

				t=b[j];
				b[j]=b[j+1];
				b[j+1]=t;

			}
			
		}
	}
}






