#ifndef _EXCEPTION_OBJ_
#define _EXCEPTION_OBJ_

/*
An object that is used to represent exceptions.
*/

class t_exception_obj
{
public:
	char* exception_msg;
	t_exception_obj(char* _exception_msg);
	~t_exception_obj();

	// Dump the exception's message.
	void panic();
};

#endif // _EXCEPTION_OBJ_