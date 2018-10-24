#include<Python.h>
#include<iostream>
using namespace std;

static PyObject *helloWorld( PyObject *self)
{
  return Py_BuildValue("s","Hello, Python extensions!!");
}

static char helloWorld_docs[] =
    "helloWorld( ): Any message you want to put here!!\n";

static PyMethodDef helloWorld_funcs[] = {
    {"helloWorld", (PyCFunction)helloWorld, 
     METH_NOARGS, helloWorld_docs},
    {NULL}
};

void inithelloWorld(void)
{
    Py_InitModule3("helloWorld", helloWorld_funcs,
                   "Extension module example!");
}
