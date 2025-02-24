/**
 * \file tokenizer.c
 * \date 3 juil. 2014
 * \author Nicolas Fedy
 */
#include <stdio.h>
#include "tokenizer.h"

int tokenizer(char* str, const char delim, char* token){

	static char* last;
	int found = 0, length = 0;
	char* check = NULL;
	int i;

	if(str)
		last = str;

	else if (!last ||*last == 0)
		return 0;

	check = last;

	while (check && !found && last[length]){

		if (*check == delim || *check == 0)
			found = 1;
		else{
			check++;
			length++;
		}
	}

	if (!found)
		return 0;

	for(i = 0; i < length; i++){
		token[i] = *last++;
	}

	token[length] = 0;
	last++;

	return 1;
}
