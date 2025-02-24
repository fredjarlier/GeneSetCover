#include "stdio.h"
#include "stdlib.h"

char in(char c, int lb, int ub)
{
	if(c > lb && c < ub)
		return 1;
	else
		return 0;
}

int main(int argc, char**argv)
{
	FILE* file = fopen(argv[1], "r");
	
	size_t len = 255;
	char* buff=malloc(len);
	char motif[len];
	motif[0]='\0';
	
	int i, j;
 	size_t index = 0;
	int offset_buff = 0;
	
	while(getline(&buff, &len, file) != -1)
	{
		offset_buff=0;
		if(motif[0]==0){
			j=2;
			motif[0] = buff[0];
			motif[1] = buff[1];
			while(buff[j]>='A' && buff[j]<='Z' && buff[j] != '\t')
			{
				motif[j]=buff[j];
				j++;
			}
			motif[j]=0;
		}
		
		else{
			int match =1;
			

			if(buff[0]!=motif[0] || buff[1]!=motif[1]){
				match = 0;
			}

			j=2;
			while(match!=0 && motif[j] != 0){
				if(buff[j]!=motif[j]){
					match = 0;
				}
				j++;
			}
			//if it doesn't match we change the motif and the index
			if(match==0){
				j=2;
				motif[0] = buff[0];
				motif[1] = buff[1];
				while(buff[j]>='A' && buff[j]<='Z' && buff[j] != '\t')
				{
					motif[j]=buff[j];
					j++;
				}
				motif[j]=0;
				index++;
			}
		}
		
		//we find the end of the line
		while(buff[offset_buff] != '\n'){
			offset_buff++;
		}
		
		buff[offset_buff] = 0;
		
		printf("%s\t%zu\n",buff,index);
	}
	free(buff);
	fclose(file);
}
