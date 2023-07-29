#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prxytypr_exception_obj.h"

t_exception_obj::t_exception_obj(char* _exception_msg)
{
	this->exception_msg = new char[strlen(_exception_msg) + 2];
	strcpy(this->exception_msg, _exception_msg);
}

t_exception_obj::~t_exception_obj()
{
	delete [] this->exception_msg;
}

void t_exception_obj::panic()
{
	printf("%s\n", this->exception_msg);
}