/* TREEDIST  */
/* Written by Chris Creevey */
/* chris.creevey@gmail.com */
/* Version 1.1  */
/* This software returns the path length distance between two taxa on the specified tree */ 

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <signal.h>
#include <time.h>

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TEMP
#define TEMP 2
#endif

time_t interval1, interval2, interval3;
double telapsed;
void unroottree(char *tree);
void weighted_pathmetric(char *string, float **scores);
int get_taxa_int(char *text);
void clear_memory(void);
void memory_error(int position);

int taxa_found =0, num_taxa=0, string_length=0;
float **scores='\0';
char *string='\0', **taxa_names ='\0';

int main(int argc, char *argv[])
    {
    int i=0, j=0, k=0, charactercount = -1, variable = 0, node_number = -1, open = 0, found = FALSE, fromroot = FALSE, first = FALSE, averagelength = FALSE, outfmt=0; 
    char number[1000], infilename[1000], outfilename[1000], c = '\0', temp[100], number2[1000];
    float *taxa_weights = '\0', *node_weights = '\0',  **closeP = '\0', branchlength = 0, totallength = 0, lastlength =0, firstlength = 0, x=0, y=0, increase = 0.1;
	FILE *infile = '\0', *outfile = '\0';
	

	if(argc < 2)
		{
		printf("treedist_all:\n\tThis sofware returns the distance between all branches of a tree.\n\tYou must call this software with the name of the file containing the tree, \n\n\tUsage: treedist_all <Tree file name> <matrix|vector> \n\n");
		printf("\tOptional output format specification: matrix (default) or vector\n");
		exit(1);
		}
	if(argc == 3)
		{
		if(strcmp(argv[2], "vector") == 0) 
			outfmt=1;
		else		
			{	
			if(strcmp(argv[2], "matrix") == 0) 
				outfmt=0;
			else
				{	
				printf("Error: output format %s unknown, please specify \"matrix\" or \"vector\"\n", argv[2]);
				exit(1);
				}
			}
		}

	number[0] = '\0';
	infilename[0] = '\0';
	temp[0] = '\0';
	strcpy(infilename, argv[1]);



	if((infile = fopen(infilename, "r")) == '\0')		/* check to see if the file is there */
		{								/* Open the source tree file */
		printf("Cannot open file %s\n", infilename);
		exit(1);
		}
	strcpy(outfilename, infilename);
	strcat(outfilename, ".dist");
	if((outfile = fopen(outfilename, "w")) == '\0')		/* check to see if the file is there */
		{								/* Open the source tree file */
		printf("Cannot open file %s\n", outfilename);
		exit(1);
		}
	
	/* read in tree to determine is length and the number of taxa */

	interval1 = time(NULL);
		
	while((c = getc(infile)) != ';' && !feof(infile))
		{
		if(c == ',') num_taxa++;
		string_length++;
		}
	num_taxa++;
	rewind(infile);
	string_length=string_length*1000;
	string=malloc((string_length)*sizeof(char));
	if(!string) memory_error(2543);
	string[0] ='\0';
	i =0;

	printf("Number of taxa in tree: %d\n", num_taxa);
	
	/* Assign the matrix to record the distances of eac htaxa to every other */
	scores=malloc(num_taxa*sizeof(float *));
	for(i=0; i<num_taxa; i++)
			{
			scores[i]=malloc(num_taxa*sizeof(float));
			for(j=0; j<num_taxa; j++)
				scores[i][j] = 0;
			}

	/*assign the array to old the names of the taxa */

	taxa_names=malloc(num_taxa*sizeof(char*));
	for (i = 0; i < num_taxa; i++)
	{
		taxa_names[i] = malloc(100*sizeof(char));
		taxa_names[i][0]= '\0';
	}
	

	/* Now read in the tree, replacing taxa names with the corresponding taxa number */
	
	i=0;
	c = getc(infile);
    while(c != ';' && !feof(infile))
        {
        switch(c)
            {
            case ')': 
				string[i] = c;
            	i++;
				while((c = getc(infile)) != '(' && c != ')' && c != ',' && c != ';' && c != ':') /* in case there are weights (like BP supports ) on the tree */
            		{
            		;
            		/*string[i] = c;
            		i++; */
            		/* don't read in the labels on the internal node, as the pathmetric algorithm assumes there are none  */
            		}
            	break;
            case ',':
            case '(':
            	string[i] = c;
            	i++;
            	c = getc(infile);
            	break;
            case ':':
            	string[i] = c;
            	i++;
            	while((c = getc(infile)) != '(' && c != ')' && c != ',' && c != ';')
            		{
            		string[i] = c;
            		i++;
            		}
            	break;
            default:
                 /*read in the taxa name, and assign it to the taxa_name array and put in the number of the taxa instead in the tree */
            	j=0;
            	number[0] = '\0';
            	number[j]=c;
            	j++;
            	while((c = getc(infile)) != '(' && c != ')' && c != ',' && c != ';' && c != ':')
            		{
					number[j]=c;
					j++;
            		}
            	number[j] = '\0';
            	number2[0] ='\0';

            	sprintf(number2, "%i", get_taxa_int(number));
            	for(j=0; j<strlen(number2); j++)
            		{
            		string[i] = number2[j];
            		i++;
            		}
                break;
            }
        }

	string[i] = ';';
	string[i+1] = '\0';
/*	printf("%s\n", string); */
	fclose(infile);

	interval2 = time(NULL);
	printf("\tFinished reading tree. Time taken %d seconds\n", (int)difftime(interval2, interval1));
	fflush(stdout);
	/* The array characters is used to keep track, for each taxa, the open and closed brackets that has followed each */
   /*	unroottree(string); */
    i=0;
    printf("Starting Pathmetric\n");
    weighted_pathmetric(string, scores); 
    fflush(stdout);
    interval3 = time(NULL);
	printf("\tFinished calculating Pathmetrics. Time taken %d seconds\n", (int)difftime(interval3, interval2));
    printf("Started printing distances\n");
    fflush(stdout);

    if(outfmt == 0)
    	{
	    for (i=0; i<num_taxa; i++)
	    	fprintf(outfile, "\t%s", taxa_names[i]);
	   	fprintf(outfile, "\n");
	    for (i=0; i<num_taxa; i++)
	    	{
	    	fprintf(outfile, "%s\t", taxa_names[i]);
   	 		for(j=0; j<num_taxa; j++)
    			{
    			fprintf(outfile, "%f\t", scores[i][j]);
    			}
    		fprintf(outfile, "\n");
    		}
    	}
    else
    	{
    	for (i=0; i<num_taxa; i++)
    		{
    		for(j=i+1; j<num_taxa; j++)
    			{
    			fprintf(outfile, "%s\t%s\t%f\n", taxa_names[i], taxa_names[j], scores[i][j]);
    			}
    		}
    	}
	
    fclose(outfile);

    interval2 = time(NULL);
	printf("\tFinished printing distance to file. Time taken %d seconds\n", (int)difftime(interval2, interval3));

    printf("Starting Clearing Memory\n");
    fflush(stdout);
    clear_memory();
    interval2 = time(NULL);
	printf("\tFinished. Overall time taken %d seconds\n", (int)difftime(interval2, interval1));
	fflush(stdout);
	}

void clear_memory(void)
	{
	int i=0;
	if(string != '\0') 
		free(string);
	if(scores!= '\0')
		{
		for(i=0; i<num_taxa; i++)	
			free(scores[i]);
		free(scores);
		}
	if(taxa_names!= '\0')
		{
		for (i = 0; i < num_taxa; i++)
			free(taxa_names[i]);
		free(taxa_names);
		}

	}

	
int get_taxa_int(char *text)
	{
	int found=FALSE, result=-1, i=0;
	/* convert string into taxa number */
	for (i = 0; i < taxa_found; i++)
		{	
		if(strcmp(text, taxa_names[i]) == 0)
			{
			found = TRUE;
			result = i;
			}
		}
	if(!found)
		{
		strcpy(taxa_names[taxa_found], text);
		result=taxa_found;
		taxa_found++;
		}
	return(result);
	}

void memory_error(int position)
	{
	interval2 = time(NULL);
	printf(" ERROR out of memory at %d, time taken up to error: %e seconds\n", position, difftime(interval2, interval1));
	clear_memory();
	exit(1);
	}


void weighted_pathmetric(char *string, float **scores)
    {
    int i=0, j=0, k=0, l=0, m=0, charactercount = -1, variable = 0, node_number = -1, open = 0, lentmpstring=0, lenstring=0, **numbrackets = '\0';
    char number[30], *tmpstring ='\0';
    float *taxa_weights = '\0', *node_weights = '\0',  **closeP = '\0';


    tmpstring=malloc(string_length*sizeof(char));
    if(!tmpstring) memory_error(444);
    tmpstring[0] ='\0';

     
    /* The array characters is used to keep track, for each taxa, the open and closed brackets that has followed each */
    
	closeP = malloc((num_taxa)*(sizeof(float*)));
    if(!closeP) memory_error(90);
    numbrackets= malloc((num_taxa)*(sizeof(int*)));
    if(!numbrackets) memory_error(908);
    for(i=0; i<num_taxa; i++)
        {
        closeP[i] = malloc(3*sizeof(float));
        if(!closeP[i]) memory_error(91);
        numbrackets[i] = malloc(2*sizeof(int));
        if(!numbrackets) memory_error(909);
        closeP[i][0] = numbrackets[i][0] = 0;  /* the number of open parentheses */
        closeP[i][1] = numbrackets[i][1] = 0;  /* the number of close parentheses */
        closeP[i][2] = FALSE;  /* whether this taxa has been found yet */
        }

    taxa_weights = malloc(num_taxa * sizeof(float));
        if(!taxa_weights) memory_error(92);
    node_weights = malloc(num_taxa * sizeof(float));
        if(!node_weights) memory_error(93);
    for(i=0; i<num_taxa; i++)
        {
        taxa_weights[i] = 0;
        node_weights[i] = 0;
        }
    unroottree(string);
   /* printf("unrooted:\n%s\n", string); */
    i=0;
    while(string[i] != ';')  /* until the end of the string */
        {
        switch(string[i])
            {
            case '(':	
                        if(i != 0)
                            {
                            k=i;
                            open = 1;
                            do
                                {
                                k++;
                                if(string[k] == '(') open++;
                                if(string[k] == ')') open--;
                                
                                }while(string[k] != ')' || open != 0);
                            /* insert into the tree after the close bracket the numberof the node for future reference */
							node_number++;
							number[0] = '\0';
							sprintf(number, "%d", node_number);
							for(l=0; l<=k; l++)
									tmpstring[l]=string[l];
							tmpstring[l] = '\0';
							

							/*string_length=string_length+strlen(number)+10;

							string = realloc(string, string_length*sizeof(char));
							if(string == '\0') memory_error(655);

							tmpstring = realloc(tmpstring, string_length*sizeof(char));
							if(tmpstring == '\0') memory_error(654); */
							
							strcat(tmpstring, number);
							l=strlen(tmpstring);
							m=k+1;
							while(string[m] != ';')
								{
								tmpstring[l] = string[m];
								l++; m++;
								}
							tmpstring[l] = ';';
							tmpstring[l+1] = '\0';
							strcpy(string, tmpstring);
							k+=strlen(number);
							/***/

							if(string[k+1] == ':')
								{
								k = k+2;  /* skip over the ':' to the start of the number */
								for(j=0; j<30; j++) number[j] = '\0';  /* initialise the array to hold the number in text form */
								j=0;
								while(string[k] != ')' && string[k] != ',' && string[k] != '(' && string[k] != ';')
									{
									number[j] = string[k];
									k++; j++;
									}  /* read in the weight */
								/*printf("%s\t%f\n", number, atof(number)); */
								node_weights[node_number] = atof(number);
								} /* otherwise the software will assume branch lengths of 1 for all */
                           
                           
                            for(j=0; j<num_taxa; j++)
                                {
                                if(closeP[j][2] == TRUE)
                                    {
                                    numbrackets[j][0] ++;
                                    closeP[j][0] += node_weights[node_number];
                                    }
                                }
                            }
                        i++;
                        break;
            case ')':	
                        if(string[i+1] != ';')
                            {
                            i++;
                            /* read in the node number */
                            l=0;
                            while(string[i] != ')' && string[i] != '(' && string[i] != ',' && string[i] != ';' && string[i] != ':')
                                {
                                number[l] = string[i];
                                i++; l++;
                            	}
                            number[l] = '\0';
                            variable=atoi(number);
                         /*   printf("variable = %d\n", variable); */
                            /* skip past the weight */
                            l=0;
                            while(string[i] != ')' && string[i] != '(' && string[i] != ',' && string[i] != ';' )
                            	{
                            	number[l] = string[i];
                            	i++; l++;
                            	}
                            number[l] = '\0';
                          /*  printf("variable = %d\ttext=%s\tweight=%f\n", variable, number, node_weights[variable]); */
                            for(j=0; j<num_taxa; j++)
                                {
                                if(closeP[j][2] == TRUE)
                                    {
                                    if(numbrackets[j][0] > 0)
                                    	{		
                                        closeP[j][0] -= node_weights[variable];  /* if this close parenthesis cancels out a previously counted openparentheis */
                                    	numbrackets[j][0] --;
                                    	}
                                    else
                                    	{	
                                        closeP[j][1] += node_weights[variable];  
                                        numbrackets[j][1]++;
                                    	}
                                    }
                                }
                            node_weights[variable] = 0; 
                                                       
                            }
                        else
                            i++;
                        break;
            case ',':
                        i++;
                        break;
            default:
                        /* this has to be a taxa number */
                        for(j=0; j<30; j++) number[j] = '\0';  /* initialise the array to hold the number in text form */ 
                        j=0;
                        while(string[i] != ':' && string[i] != ')' && string[i] != '(' && string[i] != ',')
                            {
                            number[j] = string[i];
                            i++; j++;
                            }

                        /* now need to change the string to an integer number */
                        /*printf("%s\t%d\n", number, atoi(number)); */
                        variable=atoi(number);
                        closeP[variable][2] = TRUE;  /* tell the program that this caracter has now been passed and is to be counted from now on */
                        /* Now read in the weight for this taxa */
						if(string[i] == ':')
							{
							for(j=0; j<30; j++) number[j] = '\0';  /* initialise the array to hold the number in text form */
							j=0; i++;
							while(string[i] != ')' && string[i] != ',' && string[i] != '(' && string[i] != ';')
								{
								number[j] = string[i];
								i++; j++;
								}
							/*printf("%s\t%f\n", number, atof(number));*/
							taxa_weights[variable] = atof(number);
							}
						
                        /* now need to assign the distance from any taxa we passed previously to this taxa */
                        
                        for(j=0; j<num_taxa; j++)
                            {
                            if(closeP[j][2] == TRUE && j != variable) /* we have passed this taxa earlier */
                                {
                                scores[j][variable] += (closeP[j][0]+closeP[j][1]+taxa_weights[variable]+taxa_weights[j]);  /* the score is equal to the number of open parentheses not cancelled out plus the number of close parentheses not cancelled out + 1 */
                                scores[variable][j] = scores[j][variable];
                                }
                            }
           
                        break;
            }
            
    
        }

 /*   printf("%s\n", string); */
	/*for(i=0; i<num_taxa; i++)
		printf("%f\t", taxa_weights[i]);
	printf("\n");
	for(i=0; i<num_taxa; i++)
		printf("%f\t", node_weights[i]);
	printf("\n");
	*/
/*	printf("string_length = %d\n", string_length); */

    for(i=0; i<num_taxa; i++)
    	{
        free(closeP[i]);
        free(numbrackets[i]);
        }
       free(numbrackets);
    free(closeP);
    closeP = '\0';
	free(taxa_weights);
	free(node_weights);
	free(tmpstring);

    }

  void unroottree(char * tree)
    {
    int i=0, j=0, k=0, l=0, m=0, basecount = 0, parentheses=0;
    int foundopen = FALSE, foundclose = FALSE, rooted=FALSE;
	float del_nodelen = 0;
	char length[100], *restof='\0';
	
/*	do { 
		
		rooted=FALSE; */
		restof=malloc(string_length*sizeof(char));
		i=0; j=0; k=0; l=0; m=0; basecount = 0; parentheses=0;
		foundopen = FALSE; foundclose = FALSE;
		del_nodelen = 0;
		restof[0] = '\0';
		length[0] = '\0';
		/* scan through the tree counting the number of taxa/nodes at the base (for it to be unrooted there should be at least three) */
		while(tree[i] != ';')
			{
			switch(tree[i])
				{
				case '(':
						parentheses++;
						i++;
						break;
				case ')':
						parentheses--;
						i++;
						break;
				case ',':
						if(parentheses == 1)
							{
							basecount++;
							}
						i++;
						break;
				default:
						i++;
						break;
				}
			}
			
		if(basecount <2)  /* if the base of the tree is rooted */
			{
			rooted=TRUE;
			i=0;
			parentheses = 0;
			while(tree[i] != ';')  /* delete the two parentheses to make the tree unrooted */
				{
				switch(tree[i])
					{
					case '(':
							parentheses++;
							if(parentheses == 2 && !foundopen)
								{
								tree[i] = '^';
								foundopen = TRUE;
								}
							i++;
							break;
					case ')':
							if(parentheses == 2 && !foundclose)
								{
								tree[i] = '^';
								i++;
								while(tree[i] != ')' && tree[i] != '(' && tree[i] != ',' && tree[i] != ';' && tree[i] != ':')
									{
									tree[i] = '^';
									i++;
									}
								if(tree[i] == ':')
									{
									k=0;
									length[0] = '\0';
									while(tree[i] != ')' && tree[i] != '(' && tree[i] != ',' && tree[i] != ';')
										{
										if(tree[i] != ':')
											{
											length[k] = tree[i];
											k++;
											}
										tree[i] = '^';
										i++; 
										}
									length[k] = '\0';
									}
								if(length[0] != '\0') /* we have a branch length on the internal branch, so we need to add it to the branch length of the component that is to the direct right of this closing parenthesis */
									{
									del_nodelen = atof(length);
									k=i+1; /* This should be whatever is after the ',' which should be the next compoonent */
									if(tree[k] == '(') /* we need to find the end of this clade and add the value there */
										{
										l=1; k++;
										while((l != 0 || tree[k-1] != ')') && tree[k] != ';' )   /* CHANGED RECENTLY FROM while(l != 0 && tree[k-1] != ')' && tree[k] != ';' ) */
											{
											switch(tree[k])
												{
												case '(':
													l++;
													k++;
													break;
												case ')':
													l--;
													k++;
													break;
												default:
													k++;
													break;
												}
											}
										k--; /* k now points to the closing bracket */
										/* read in the length attached to this partenthsis */
										k++;
										while(tree[k] != ')' && tree[k] != '(' && tree[k] != ',' && tree[k] != ';' && tree[k] != ':') k++; /* k now points to the ":" at the start of the length (if there is a length) */
										}
									else
										{
										while(tree[k] != ')' && tree[k] != '(' && tree[k] != ',' && tree[k] != ';' && tree[k] != ':') k++; /* k now points to the ":" at the start of the length (if there is a length) */
										}
									if(tree[k] == ':') /* there is length attached to this */
										{
										m=k+1;
										length[0] = '\0'; l=0;
										while(tree[m] != ')' && tree[m] != '(' && tree[m] != ',' && tree[m] != ';')
											{
											length[l] = tree[m];
											l++; m++;
											}
										length[l] = '\0';
										del_nodelen += atof(length);
										
										}
									else
										m=k;
									/* now add del_nodelen to this point in the tree */
									 l=0;
									while(tree[m] != ';' && tree[m] != '\0')
										{
										restof[l] = tree[m];
										m++; l++;
										}
									restof[l] = ';';
									restof[l+1] = '\0';
									if(tree[k] == ':')
										tree[k] = '\0';
									else
										{
										tree[k] = '\0';
										}
									length[0] = '\0';
									sprintf(length, ":%f", del_nodelen);
									strcat(tree, length);
									strcat(tree, restof);
									}
								foundclose = TRUE;
								}
							i++;
							parentheses--;
							break;
					default:
							i++;
							break;
					}
				}
							
			/* scan through the string shifting up the characters to take into account those parentheses that have been deleted */
			i=0; j=0;
			while(tree[j] != ';')
				{
				if(tree[j] == '^')
					{
					while(tree[j] == '^')
						j++;
					if(i!= j)tree[i] = tree[j];
					i++; j++;
					}
				else
					{
					if(i!=j)tree[i] = tree[j];
					i++;j++;
					}
				}
			tree[i] = tree[j];
			tree[i+1] = '\0';
			
			}
	/*	} while(rooted==TRUE); */
	free(restof);
    } 

