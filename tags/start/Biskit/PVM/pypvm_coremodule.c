#include <Python.h>
#include <pvm3.h>
#include <string.h>
#include <stdlib.h>

/*    FILE: pypvm_core_module.c -- A python interface to PVM.  For more
 *         information on PVM, see http://www.epm.ornl.gov/pvm/.
 * AUTHORS: W. Michael Petullo, wp0002@drake.edu
 *          Gregory D. Baker,  greg.baker@ifost.org.au
 *    DATE: 13 FEB 1998 - 24 MAY 2000
 *
 * Copyright (c) 1999 W. Michael Petullo
 * Copyright (c) 2000 Gregory D. Baker
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* Useful definitions for this module */

static PyObject * pypvm_module;
static PyObject * pypvm_dictionary;

static PyObject 
  *PvmBadParamException, *PvmMismatchException,	 *PvmOverflowException,	 
  *PvmNoDataException ,	 *PvmNoHostException,	 *PvmNoFileException,	 
  *PvmNoMemException ,	 *PvmBadMsgException,	 *PvmSysErrException,	 
  *PvmNoBufException ,	 *PvmNoSuchBufException, *PvmNullGroupException,
  *PvmDupGroupException, *PvmNoGroupException,   *PvmNotInGroupException ,
  *PvmNoInstException ,	 *PvmHostFailException,  *PvmNoParentException ,	 
  *PvmNotImplException , *PvmDSysErrException,	 *PvmBadVersionException,	 
  *PvmOutOfResException, *PvmDupHostException,	 *PvmCantStartException,	 
  *PvmAlreadyException , *PvmNoTaskException,	 *PvmNoEntryException,	 
  *PvmDupEntryException,


  /* one more for good luck  ;-) */
  *pypvmUnknownExceptionException
  ;

static PyObject 
  *PvmBadParamNumber, *PvmMismatchNumber,	 *PvmOverflowNumber,	 
  *PvmNoDataNumber ,	 *PvmNoHostNumber,	 *PvmNoFileNumber,	 
  *PvmNoMemNumber ,	 *PvmBadMsgNumber,	 *PvmSysErrNumber,	 
  *PvmNoBufNumber ,	 *PvmNoSuchBufNumber, *PvmNullGroupNumber,
  *PvmDupGroupNumber, *PvmNoGroupNumber,   *PvmNotInGroupNumber ,
  *PvmNoInstNumber ,	 *PvmHostFailNumber,  *PvmNoParentNumber ,	 
  *PvmNotImplNumber , *PvmDSysErrNumber,	 *PvmBadVersionNumber,	 
  *PvmOutOfResNumber, *PvmDupHostNumber,	 *PvmCantStartNumber,	 
  *PvmAlreadyNumber , *PvmNoTaskNumber,	 *PvmNoEntryNumber,	 
  *PvmDupEntryNumber;


/* There are lots of functions that may want to return nothing,  or
   throw an exception if their corresponding PVM call returned < 0.
   We define a macro `return_none' to handle this.  Its argument is
   the return result from the PVM call.  (Which will probably be
   a variable name. */

#define return_none(x) \
  if (was_error(x)) return NULL; \
  Py_INCREF(Py_None); \
  return Py_None;

/* PVM often returns negative numbers to mean some kind of error.
   Python can throw exceptions in these cases,  so that's what
   we do.  (We use set object with the exception as the string,
   and the data as the PVM return value.  Maybe it should be
   the other way around.  Anyway,  it's kind of cute. */


static int was_error(int info) {
  if (info > 0) return 0;
  switch (info) {
  case PvmOk: return 0;
#define error_case(x) case x : \
    PyErr_SetObject(x ## Exception,x ## Number); \
    return 1;
  error_case(PvmBadParam); error_case(PvmMismatch);  error_case(PvmOverflow);
  error_case(PvmNoData);  error_case(PvmNoHost);  error_case(PvmNoFile);
  error_case(PvmNoMem);  error_case(PvmBadMsg);  error_case(PvmSysErr);
  error_case(PvmNoBuf); error_case(PvmNoSuchBuf);  error_case(PvmNullGroup);
  error_case(PvmDupGroup);  error_case(PvmNoGroup);  error_case(PvmNotInGroup);
  error_case(PvmNoInst);  error_case(PvmHostFail);  error_case(PvmNoParent);
  error_case(PvmNotImpl);  error_case(PvmDSysErr);  error_case(PvmBadVersion);
  error_case(PvmOutOfRes);  error_case(PvmDupHost);  error_case(PvmCantStart);
  error_case(PvmAlready);  error_case(PvmNoTask);  error_case(PvmNoEntry);
  error_case(PvmDupEntry);
  }

  PyErr_SetObject(pypvmUnknownExceptionException,
		  PyInt_FromLong( (long) info));
  return 1;
}


static char pypvm_hostinfo__doc[] =
"Returns a list (one tuple for each host in the virtual machine) - each tuple\n"
"is (task id of pvm daemon,hostname,architecture,speed).\n"
"This is a subset of the information that can be retrieved with pypvm.config()\n"
;

static PyObject * pypvm_hostinfo (PyObject *self, PyObject *args, PyObject* keywords) {
  struct pvmhostinfo *hinfo;
  int i;
  PyObject *list;
  int nhost;
  int narch;
  int info;

  PyObject *pvmd_tid, *hostname, *architecture, *speed;
  PyObject *tuple;

  info = pvm_config(&nhost,&narch,&hinfo);
  if (was_error(info)) return NULL;
  list = PyList_New(nhost);

  for (i = 0; i < nhost; i++) {
    pvmd_tid = PyInt_FromLong(hinfo[i].hi_tid);
    hostname = PyString_FromString(hinfo[i].hi_name);
    architecture = PyString_FromString(hinfo[i].hi_arch);
    speed = PyInt_FromLong(hinfo[i].hi_speed);
    
    tuple = PyTuple_New(4);
    PyTuple_SetItem(tuple,0,pvmd_tid);
    PyTuple_SetItem(tuple,1,hostname);
    PyTuple_SetItem(tuple,2,architecture); 
    PyTuple_SetItem(tuple,3,speed);
    PyList_SetItem(list,i,tuple);
  }
  return list;
}



static char pypvm_narch__doc[]= 
"pypvm.narch() returns a count of the number of different data formats\n"
"being used (narch = number of architectures).  This is a subset of\n"
"the information returned by pypvm.config()\n";

static PyObject * pypvm_narch(PyObject * self, PyObject * args, PyObject * keywords) {
  struct pvmhostinfo *hinfo;
  PyObject *list;
  int nhost;
  int narch;
  int info;

  info = pvm_config(&nhost,&narch,&hinfo);
  if (was_error(info)) return NULL;
  return Py_BuildValue("i",nhost);
}




static char pypvm_config__doc[] = 
"pypvm.config() returns a three element tuple\n"
" - the number of hosts in the PVM cluster\n"
" - the number of architectures in the cluster\n"
" - a list of dictionaries giving information about each machine in\n"
"   the cluster.\n"
"\n"
" Each dictionary has keys \n"
" - dtid\n"
" - hostname\n"
" - arch\n"
" - speed\n"
"\n"
"Note that the speed field may have no particularly relationship\n"
"with reality,  as it just reporting the number stored in the\n"
"startup hostsfile for the master pvmd.  Read the pvmd3 man\n"
"pages for more information.\n"
"\n"
"The information provided by pypvm.config is also available throug\n"
"pypvm.hostinfo and pypvm.narch.";

static PyObject * pypvm_config (PyObject *self, PyObject *args, PyObject * keywords) {
  int info;
  int nhosts, narchs, i;
  struct pvmhostinfo *hinfo;
  PyObject *py_info, *items, *fctval;
  info =pvm_config(&nhosts, &narchs, &hinfo);
  if (was_error(info)) { return NULL; }
  py_info = PyList_New(nhosts);
  for (i = 0; i < nhosts; i++) {
    /* build info dictionaries */
    items = PyDict_New();
    PyDict_SetItemString(items, "dtid",PyInt_FromLong(hinfo[i].hi_tid));
    PyDict_SetItemString(items, "hostname",PyString_FromString(hinfo[i].hi_name));
    PyDict_SetItemString(items, "arch",PyString_FromString(hinfo[i].hi_arch));
    PyDict_SetItemString(items, "speed",PyInt_FromLong(hinfo[i].hi_speed));
    PyList_SetItem(py_info, i, items);
  }
  fctval = PyTuple_New(3);
  PyTuple_SetItem(fctval, 0, PyInt_FromLong(nhosts));
  PyTuple_SetItem(fctval, 1, PyInt_FromLong(narchs));
  PyTuple_SetItem(fctval, 2, py_info);
  free(hinfo);
  return (fctval);
}


static char pypvm_tasks__doc[] =
"Returns a list (one tuple for each task) of all tasks running in the \n"
"virtual machine.  Each tuple is (task id, parent task id, task id of\n"
"the pvm daemon on the tasks machine,flags it was called with,\n"
"the name of the task, and the process id of the task).\n"
"\n"
"If called with a tid as an argument,  it will just return that tasks\n"
"info;  if called with the tid of a pvm daemon,  it will return all\n"
"the tasks on that host.\n";

static PyObject * pypvm_tasks (PyObject* self,PyObject* args, PyObject * keywords) {
  struct pvmtaskinfo * taskinfo;
  int where=0;   /* If they don't provide it, they probably meant everything */
  int ntasks;
  int info;
  int i;

  PyObject * list, *tuple;
  static char *kwlist[] = { "tid" , NULL };

  
  if (!PyArg_ParseTupleAndKeywords(args,keywords,"|i", kwlist,&where))  return NULL;

  info = pvm_tasks(where,&ntasks,&taskinfo);
  if (was_error(info)) return NULL;
  
  list = PyList_New(ntasks);
  for (i=0;i<ntasks;i++) {
    tuple = Py_BuildValue("(iiiisi)",
			  taskinfo[i].ti_tid,
			  taskinfo[i].ti_ptid,
			  taskinfo[i].ti_host,
			  taskinfo[i].ti_flag,
			  taskinfo[i].ti_a_out,
			  taskinfo[i].ti_pid);
    PyList_SetItem(list,i,tuple);
  }

  return list;
}


/* -------------------------------------------------------------------------- */
/* The following functions were originally created by swig, but I've modified them all */
/* -------------------------------------------------------------------------- */


static char pypvm_addhosts__doc[] =
"\n"
"pypvm.addhosts(hostlist) adds the computers named in [hostlist] to the\n"
"configuration of computers making up the virtual machine.  The names\n"
"should have the same syntax as lines of a pvmd hostfile (see man page\n"
"for pvmd3): A hostname followed by options of the form xx=y.\n"
"\n"
"pypvm.addhosts returns a list;  this list will (hopefully) just contain\n"
"zeros.  However, if any of them are negative numbers,  they are \n"
"indicative that the addhost operation failed on that corresponding\n"
"machine.  The return results are equal to members of of the pypvm.results\n"
"dictionary.\n"
"  pypvm.results[BadParam]  =>     bad hostname syntax.\n"
"  pypvm.results[NoHost]    =>     no such host.\n"
"  pypvm.results[CantStart] =>     failed to start pvmd on host.\n"
"  pypvm.results[DupHost]   =>     host already configured.\n"
"  pypvm.results[BadVersion]=>     pvmd protocol versions dont match.\n"
"  pypvm.results[OutOfRes]  =>     PVM has run out of system resources.\n"
;

static PyObject * pypvm_addhosts(PyObject *self, PyObject *args, PyObject *keywords) {
    int num_hosts, i, *infos;
    char **hosts = (char **) NULL;
    PyObject *host_list, *host, *info_list;
    static char * kwlist[] = {"hostlist", NULL};
    int return_info;

    if (!PyArg_ParseTupleAndKeywords(args, keywords,"O:hostlist", kwlist, &host_list))
	/* ASSERT: incorrect arguments */
	return (NULL);
    if ((num_hosts = PyList_Size(host_list)) < 0) {
	/* ASSERT: second argument was not a list */
	PyErr_SetString(PyExc_TypeError,
			"argument 1: expected list of strings");
	return (NULL);
    }
    if (num_hosts)
	/* allocate memory for host array */
	if ((hosts = (char **) PyMem_Malloc(num_hosts * sizeof(char *) + 1))
	    == NULL)
	    return (NULL);
    for (i = 0; i < num_hosts; i++) {
	/* build host array */
	host = PyList_GetItem(host_list, i);
	if ((hosts[i] = PyString_AsString(host)) == NULL) {
	    PyErr_SetString(PyExc_TypeError,
			    "argument 1: expected list of strings");
	    return (NULL);
	}
    }
    if ((infos = (int *) PyMem_Malloc(num_hosts * sizeof(int))) == NULL)
	return (NULL);
    return_info = pvm_addhosts(hosts, num_hosts, infos);
    if (was_error(return_info)) { return NULL; }

    info_list = PyList_New(num_hosts);
    for (i = 0; i < num_hosts; i++)
	PyList_SetItem(info_list, i, PyInt_FromLong(infos[i]));
    PyMem_Free(hosts);
    PyMem_Free(infos);
    return (info_list);
}


static char pypvm_archcode__doc[] =
"The routine pypvm.archcode() returns an integer given an architecture\n"
"name.  The code returned identifies machines with compatible binary\n"
"data formats.  For example, SUN4 and RS6K have the same code, while\n"
"ALPHA has a different one (because a few datatypes have different\n"
"sizes).  This lets you know when you can get away with using\n"
"pypvm.data.raw instead of the default encoding to pass messages between\n"
"tasks on two machines.\n"
"\n"
"Naturally, you shouldnt assume the values returned by pvm_archcode\n"
"are etched in stone; the numbers have no intrinsic meaning except that\n"
"if two different arch names map to the same value then theyre\n"
"compatible.  ";

static PyObject *pypvm_archcode(PyObject *self,PyObject *args, PyObject * keywords) 
{ 
  PyObject * resultobj; 
  int result; 
  char * arch;

  static char * kwlist[] = { "arch" , NULL };

  self = self;
  if(!PyArg_ParseTupleAndKeywords(args,keywords,"s:pvm_archcode",kwlist, &arch)) 
    return NULL;
  result = (int )pvm_archcode(arch);
  resultobj = Py_BuildValue("i",result);
  return resultobj;
}


static char pypvm_barrier__doc[] =
"\n"
"pypvm.barrier(group,count)\n"
"\n"
"The routine pypvm.barrier blocks the calling process until [count] members\n"
"of the [group] have called pypvm.barrier.\n"
;


static PyObject *pypvm_barrier(PyObject *self, PyObject *args, PyObject * keywords) {
    int  info;
    char * group;
    int  count;

    static char * kwlist[] = {"group", "count" , NULL };

    self = self;
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"si:pvm_barrier",kwlist, &group,&count)) 
        return NULL;
    info = pvm_barrier(group,count);
    return_none(info);
}

static char pypvm_bcast__doc[] = 
"pypvm.bcast(group,msgtag)\n"
"\n"
"Broadcasts the data in the active message buffer to the processes in [group] \n"
;

static PyObject *pypvm_bcast(PyObject *self, PyObject *args, PyObject * keywords) {
    int  info;
    char * group;
    int  msgtag;

    static char * kwlist[] = {"group", "msgtag", NULL };

    self = self;
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"si:pvm_bcast",kwlist, &group,&msgtag)) 
        return NULL;
    info = pvm_bcast(group,msgtag);
    return_none(info);
}

static char pypvm_bufinfo__doc[] =
"pypvm.bufinfo(bufid) returns a tuple of information about a message buffer -\n"
"(bytes,msgtag,tid) where \n"
"\n"
" [bytes] is an integer returning the length in bytes of the entire message.\n"
" [msgtag] is an integer returning the message label.  Useful when the message\n"
"was received with a wildcard msgtag\n"
" [tid] is an integer returning the source of the message.  Useful when the \n"
"message was received with a wildcard tid.\n"
;

static PyObject *pypvm_bufinfo(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  info;
    int  bufid;
    int  bytes, msgtag, source_tid;
    static char * kwlist[] = {"bifid", NULL};

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"i:pvm_bufinfo",kwlist, &bufid)) 
        return NULL;
    info = pvm_bufinfo(bufid,&bytes,&msgtag,&source_tid);
    if (was_error(info)) return NULL;
    resultobj = Py_BuildValue("(iii)",bytes,msgtag,source_tid);

    return resultobj;
}


/* The pvm3.h header file seems to indicate that this function is only
   available on systems for which EOF is defined. */

static char pypvm_catchout__doc[] = 
"pypvm.catchout(file) causes the calling task (the parent)  to catch output \n"
"from tasks spawned after the call to pypvm.catchout().  Characters printed \n"
"on stdout or stderr in children tasks are collected by the pvmds and sent \n"
"in control messages to the parent task, which tags each line and appends  \n"
"it to the specified file.  Output from grandchildren (spawned by children) \n"
"tasks is also  collected,  provided the children dont reset PvmOutputTid.\n"
"\n"
"If [file] is not specified,  or is None,  output collection is turned off.\n"
"\n"
"If pypvm.exit() is called while output collection is in effect, it will block \n"
"until all tasks sending it output have exited, in order to print all their \n"
"output.  To avoid this, turn collection off first.";


/* ' */
static PyObject *pypvm_catchout(PyObject * self, PyObject * args, PyObject *keywords)
{
  int info;
  FILE *file;
  PyObject *py_file=NULL;
  static char * kwlist[] = {"file", NULL};
  if (!PyArg_ParseTupleAndKeywords(args,keywords, "|O", kwlist, &py_file))
    return (NULL);
  if ((py_file == NULL) || (py_file == Py_None)) {
    /* Obviously they want collection turned off */
    file = NULL;
  } else if (!PyFile_Check(py_file)) {
    PyErr_SetString(PyExc_TypeError, "argument 1: expected file");
    return (NULL);
  } else {
    file = PyFile_AsFile(py_file);
  }
  info = pvm_catchout(file);
  return_none(info);
}



static char pypvm_delete__doc[]=
"\n"
"pypvm.delete(name,index)  deletes entry <name, index> from the database ;\n"
"see pvm_insert(3PVM) for a description of this database. \n"
;

static PyObject *pypvm_delete(PyObject *self, PyObject *args, PyObject * keywords) {
    int  cc;
    char * name;
    int  index;
    static char * kwlist[] = {"name","index",NULL};

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"si:pvm_delete",kwlist,&name,&index)) 
        return NULL;
    cc = pvm_delete(name,index);
    return_none(cc);
}

static char pypvm_delhosts__doc[] = 
"pypvm.delhosts(hostlist) deletes the computers pointed to in [hostlist]\n"
"from the existing configuration of computers mak­ ing up the virtual\n"
"machine.  All PVM processes and the pvmd running on these computers\n"
"are killed as the computer is deleted.\n"
"\n"
"The returned array can checked to determine which host caused the error,\n"
"and compared against the elements in the pypvm.results dictionary.\n"
;

static PyObject * pypvm_delhosts(PyObject *self, PyObject* args, PyObject *keywords) {
  int return_result;
  int num_hosts, i, *infos;
  char **hosts = (char **) NULL;
  PyObject *host_list, *host, *info_list;
  char * kwlist[] = {"hostlist" , NULL };
  if (!PyArg_ParseTupleAndKeywords(args, keywords, "O", kwlist,  &host_list))
	/* ASSERT: incorrect arguments */
	return (NULL);

  if ((num_hosts = PyList_Size(host_list)) < 0) {
    /* ASSERT: second argument was not a list */
    PyErr_SetString(PyExc_TypeError,"argument 1: expected list of strings");
    return (NULL);
  }
  if (num_hosts)
    /* allocate memory for host array */
    if ((hosts = (char **) PyMem_Malloc(num_hosts * sizeof(char *) + 1))
	== NULL)
      return (NULL);
  for (i = 0; i < num_hosts; i++) {
    /* build host array */
    host = PyList_GetItem(host_list, i);
    if ((hosts[i] = PyString_AsString(host)) == NULL) {
      PyErr_SetString(PyExc_TypeError,
		      "argument 1: expected list of strings");
      return (NULL);
    }
  }
  if ((infos = (int *) PyMem_Malloc(num_hosts * sizeof(int))) == NULL)
    return (NULL);
  return_result = pvm_delhosts(hosts, num_hosts, infos);
  if (was_error(return_result)) { return NULL; }
  info_list = PyList_New(num_hosts);
  for (i = 0; i < num_hosts; i++)
    PyList_SetItem(info_list, i, PyInt_FromLong(infos[i]));
  PyMem_Free(hosts);
  PyMem_Free(infos);
  return (info_list);
}

static char pypvm_exit__doc[]=

"The routine pypvm.exit() tells the local pvmd that this process is leaving PVM.  \n"
"This routine does not kill the process, which can continue to perform tasks just \n"
"like any other serial process.\n"
;

static PyObject *pypvm_exit(PyObject *self, PyObject *args, PyObject * keywords) {
    int  result;

    result = pvm_exit();
    return_none(result);
}

static char pypvm_export__doc[]=

"pypvm.export(name) is provided for convenience in editing the environment  \n"
"variable PVM_EXPORT, while maintaining the colon-separated list syntax it\n"
"requires.\n"
;

static PyObject *pypvm_export(PyObject *self, PyObject *args, PyObject * keywords) {
  int  cc;
  char * name;
  static char * kwlist[] = {"name",NULL};

  if(!PyArg_ParseTupleAndKeywords(args,keywords,"s:pvm_export",kwlist, &name)) 
    return NULL;
  cc = pvm_export(name);
  return_none(cc);
}

static char pypvm_freebuf__doc[]= 
"The routine pypvm.freebuf( bufid ) frees\n"
"the memory associated with the message buffer identified by bufid.\n"
"Message buffers are created by pypvm.mkbuf, pypvm.initsend, and pypvm.recv.\n"
;

static PyObject *pypvm_freebuf(PyObject *self, PyObject *args, PyObject * keywords) {
    int  info;
    int  bufid;
    static char * kwlist[] = {"bufid",NULL};
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"i:pvm_freebuf", kwlist, &bufid)) 
        return NULL;
    info = pvm_freebuf(bufid);
    return_none(info);
}

/*
PyObject * pvm_gather(PyObject * self, PyObject * args) {
}
*/

/* This next function should return Python file objects.
   Better still,  we should have a python select wrapper as well. */

/*
static PyObject *pypvm_getfds(PyObject *self, PyObject *args) {
    PyObject * resultobj;
    int  result;
    int ** arg0;
    char * argc0 = 0;

    self = self;
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"s:pvm_getfds", kwlist,&argc0)) 
        return NULL;
    if (argc0) {
        if (SWIG_GetPtr(argc0,(void **) &arg0,"_int_pp")) {
            PyErr_SetString(PyExc_TypeError,"Type error in argument 1 of pvm_getfds. Expected _int_pp.");
        return NULL;
        }
    }
    result = (int )pvm_getfds(arg0);
    resultobj = Py_BuildValue("i",result);
    return resultobj;
}
*/

static char pypvm_getinst__doc[]=
"pypvm.getinst(group,tid) returns the instance number in [group] of a\n"
"PVM process identified by [tid]";


static PyObject *pypvm_getinst(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  inum;
    char * group;
    int  tid;
    static char *kwlist[] = { "group", "tid" , NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"si:pvm_getinst", kwlist,&group,&tid)) 
        return NULL;
    inum= pvm_getinst(group,tid);
    if (was_error(inum)) return NULL;
    resultobj = Py_BuildValue("i",inum);
    return resultobj;
}


#ifdef PVM33COMPAT

char pypvm_getmwid__doc[] = 
"\n"
"pypvm.getmwid(bufid) gets the wait ID of a message.\n"
"A wait identifier is part of a message (like  the  source, destination, tag and body). \n"
"It is used to match a reply to the corresponding request. The default wait ID for a \n"
"message is zero (none). \n"
;

static PyObject *pypvm_getmwid(PyObject *self, PyObject *args, PyObject * keywords) {
    int  info;
    int  bufid;
    static char *kwlist[] = { "bufid" , NULL };

    self = self;
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"i:pvm_getmwid", kwlist,&bufid))  return NULL;
    info = pvm_getmwid(bufid);
    return_none(info);
}

#endif

static char pypvm_getopt__doc[] =
"pypvm.getopt(what) returns the value of the specified option in PVM.\n"
"For a discussion of options and values, see the pvm_setopt man page.\n"
"\n"
"The permissible values of [what] are in the pypvm.opt dictionary:\n"
"  pypvm.opt[Route]            1    Message routing policy\n"
"  pypvm.opt[DebugMask]        2    Libpvm debug mask\n"
"  pypvm.opt[AutoErr]          3    Auto error reporting\n"
"  pypvm.opt[OutputTid]        4    Stdout destination for children\n"
"  pypvm.opt[OutputCode]       5    Output message tag\n"
"  pypvm.opt[TraceTid]         6    Trace data destination for children\n"
"  pypvm.opt[TraceCode]        7    Trace message tag\n"
"  pypvm.opt[FragSize]         8    Message fragment size\n"
"  pypvm.opt[ResvTids]         9    Allow messages to reserved tags and TIDs\n"
"  pypvm.opt[SelfOutputTid]   10    Stdout destination\n"
"  pypvm.opt[SelfOutputCode]  11    Output message tag\n"
"  pypvm.opt[SelfTraceTid]    12    Trace data destination\n"
"  pypvm.opt[SelfTraceCode]   13    Trace message tag\n"
"  pypvm.opt[ShowTids]        14    pvm_catchout prints task ids with output\n"
"  pypvm.opt[PollType]        15    Message wait policy (shared memory)\n"
"  pypvm.opt[PollTime]        16    Message spinwait duration\n"
"\n"
"Example:\n"
" pypvm.getopt(pypvm.opt[FragSize])\n"
;

static PyObject *pypvm_getopt(PyObject *self, PyObject *args, PyObject * keywords) {
    int what, val;
    char * kwlist [] = {"what", NULL };
    if (!PyArg_ParseTupleAndKeywords(args,keywords, "i", kwlist, &what))
	/* incorrect arguments */
	return (NULL);
    val = pvm_getopt(what);
    if (was_error(val)) { return NULL; }
    return (PyLong_FromLong(val));
}



static char pypvm_getrbuf__doc[] = 
"\n"
"pypvm.getrbuf() returns  the message buffer identifier for the active receive buffer.\n"
;

static PyObject *pypvm_getrbuf(PyObject *self, PyObject *args, PyObject * keywords) 
{
    PyObject * resultobj;
    int  result;

    if(!PyArg_ParseTuple(args,":pvm_getrbuf"))  return NULL;
    result = (int )pvm_getrbuf();
    resultobj = Py_BuildValue("i",result);
    return resultobj;
}

static char pypvm_getsbuf__doc[] =
"pypvm.getsbuf() returns  the message buffer identifier for the active send buffer.";

static PyObject *pypvm_getsbuf(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  result;

    if(!PyArg_ParseTuple(args,":pvm_getsbuf"))  return NULL;
    result = (int )pvm_getsbuf();
    resultobj = Py_BuildValue("i",result);
    return resultobj;
}


static char pypvm_gettid__doc[]=
"pypvm.gettid(group,instance) returns the tid of the process identified by\n"
"[group] and [instance].   \n"
;
static PyObject *pypvm_gettid(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  tid;
    char * group;
    int  inum;
    static char *kwlist[] = { "group","instance", NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"si:pvm_gettid", kwlist,&group,&inum)) return NULL;
    tid = pvm_gettid(group,inum);
    if (was_error(tid)) return NULL;

    resultobj = Py_BuildValue("i",tid);
    return resultobj;
}

static char pypvm_gsize__doc[]=
"\n"
"pypvm.gsize(group) returns the number of members presently in [group]. \n"
;
static PyObject *pypvm_gsize(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  result;
    char * group;
    static char *kwlist[] = {"group" , NULL };

    self = self;
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"s:pvm_gsize", kwlist,&group)) return NULL;
    result = (int )pvm_gsize(group);
    if (was_error(result)) return NULL;
    resultobj = Py_BuildValue("i",result);
    return resultobj;
}

static char pypvm_halt__doc[]="pypvm.halt() shuts down the entire PVM system.";
static PyObject *pypvm_halt(PyObject *self, PyObject *args, PyObject * keywords) {
    int  result;

    result = pvm_halt();
    return_none(result);
}

static char pypvm_initsend__doc[]=
"\n"
"pypvm.initsend(encoding) / pypvm.initsend() clears the default send buffer and \n"
"specify message encoding to be [encoding]/default.  Valid encoding values are \n"
"in pypvm.data and can be either pypvm.data[default],  pypvm.data[raw] or \n"
"pypvm.data[inplace].  I dont recommend using inplace with this version of PyPVM.";

static PyObject *pypvm_initsend(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  bufid;
    int  encoding;
    static char *kwlist[] = { "encoding" , NULL };

    encoding = PvmDataDefault;  /* if they don't provide an argument, then use XDR */
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"|i", kwlist,&encoding)) return NULL;

    bufid = pvm_initsend(encoding);
    if (was_error(bufid)) return NULL;
    resultobj = Py_BuildValue("i",bufid);
    return resultobj;
}

static char pypvm_insert__doc[]=
"\n"
"pypvm.insert(name,index,data) stores data in the simple database maintained by\n"
"the master pvmd , which can  be used to store values such as tids and make them accessible\n"
"anywhere within a virtual machine. This is useful when building an application\n"
"such as the group server, which must advertise its task id so clients can register\n"
"send messages to register.  [data] (an integer) gets stoed as the entry <[name], [index]>.\n"
;
static PyObject *pypvm_insert(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  cc;
    char * name;
    int  index;
    int  data;
    static char * kwlist[] = {"name", "index", "data", NULL};

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"sii:pvm_insert", kwlist,&name,&index,&data)) return NULL;
    cc = pvm_insert(name,index,data);
    if (was_error(cc)) return NULL;
    resultobj = Py_BuildValue("i",cc);
    return resultobj;
}


static char pypvm_joingroup__doc[]=
"\n"
"pypvm.joingroup(group) enrolls the calling process in [group]\n"
;

static PyObject *pypvm_joingroup(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  inum;
    char * group;
    static char *kwlist[] = { "group" , NULL };

    self = self;
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"s:pvm_joingroup", kwlist,&group)) return NULL;
    inum= pvm_joingroup(group);
    if (was_error(inum)) return NULL;
    resultobj = Py_BuildValue("i",inum);
    return resultobj;
}

static char pypvm_kill__doc[]=
"pypvm.kill(tid) sends a terminate (SIGTERM) signal to the PVM process [tid]. \n"
;
static PyObject *pypvm_kill(PyObject *self, PyObject *args, PyObject * keywords) {
    int  info;
    int  tid;
    static char *kwlist[] = { "tid" , NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"i:pvm_kill", kwlist,&tid)) return NULL;
    info = pvm_kill(tid);
    return_none(info);
}

static char pypvm_lookup__doc[]=
"pypvm.lookup(name,index) retrieves data stored in the location given by\n"
"<[name],[index]>.  If index is -1, the  data  stored  at  the first \n"
"existing index in the named class is returned. \n"
"\n"
"See pvm_insert(3PVM) (or pypvm.insert.__doc__) for a description of this database.  \n"
;

static PyObject *pypvm_lookup(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  cc;
    char * name;
    int  index;
    int data;
    static char *kwlist[] = { "name", "index" , NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"si:pvm_lookup", kwlist,&name,&index)) return NULL;
    cc = pvm_lookup(name,index,&data);
    if (was_error(cc)) return NULL;
    resultobj = Py_BuildValue("i",cc);
    return resultobj;
}

static char pypvm_lvgroup__doc[]=
"\n"
"pypvm.lvgroup(group) unenrolls the calling process from a named [group].  \n"
;

static PyObject *pypvm_lvgroup(PyObject *self, PyObject *args, PyObject * keywords) {
    int  info;
    char * group;
    static char *kwlist[] = {"group" , NULL };

    self = self;
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"s:pvm_lvgroup", kwlist,&group)) 
        return NULL;
    info = pvm_lvgroup(group);
    return_none(info);
}

static char pypvm_mcast__doc[] =
"pypvm.mcast(tids,msgtag) multicasts a message stored in the active send\n"
"buffer to the tasks specified in the [tids] list.  The message is not\n"
"sent to the caller even if listed in [tids].  The content of the\n"
"message can be distinguished by [msgtag].\n"
;

static PyObject *pypvm_mcast(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject *tids_list, *tid;
    int *tids, msgtag, list_len, i;
    int info;
    char * kwlist[] = {"tids", "msgtag", NULL };
    if (!PyArg_ParseTupleAndKeywords(args, keywords, "Oi", kwlist, &tids_list, &msgtag))
	return (NULL);
    if ((list_len = PyList_Size(tids_list)) < 0) {
	/* ASSERT: second argument was not a list */
	PyErr_SetString(PyExc_TypeError,
			"argument 1: expected list of tids");
	return (NULL);
    }
    if ((tids = (int *) PyMem_Malloc(list_len * sizeof(int))) == NULL)
	return (NULL);
    for (i = 0; i < list_len; i++) {
	tid = PyList_GetItem(tids_list, i);
	/* extract tids from list */
	if (!PyInt_Check(tid)) {
	    PyErr_SetString(PyExc_TypeError,
			    "argument 1: expected list of tids");
	    return (NULL);
	}
	tids[i] = PyInt_AsLong(tid);

    }
    info = pvm_mcast(tids, list_len, msgtag);
    PyMem_Free(tids);
    return_none(info);
}


static char pypvm_mkbuf__doc[]=
"\n"
"pypvm.mkbuf(encoding) creates a new message buffer. Valid encoding values are \n"
"in pypvm.data and can be either default,  raw or inplace.  I dont recommend \n"
"using inplace with this version of PyPVM.  The encoding argument is optional.\n"
;
static PyObject *pypvm_mkbuf(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  bufid;
    int  encoding;
    static char *kwlist[] = { "encoding", NULL };

    encoding = PvmDataDefault;   /* if they don't provide an encoding,  they mean XDR */
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"|i", kwlist,&encoding)) return NULL;
    bufid= pvm_mkbuf(encoding);
    if (was_error(bufid)) return NULL;
    resultobj = Py_BuildValue("i",bufid);
return resultobj;
}

static char pypvm_mstat__doc[]=
"pypvm.mstat(host) returns the status of a host in the virtual machine.\n"
"   value               MEANING\n"
"   PvmOk               host is OK\n"
"   PvmNoHost           host is not in virtual machine\n"
"   PvmHostFail         host is unreachable (and thus possibly failed)";

static PyObject *pypvm_mstat(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  mstat;
    char * host;
    static char *kwlist[] = { "host" , NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"s:pvm_mstat", kwlist,&host)) return NULL;
    mstat = pvm_mstat(host);
    if (was_error(mstat)) return NULL;
    resultobj = Py_BuildValue("i",mstat);
    return resultobj;
}

static char pypvm_mytid__doc[]="pypvm.mytid() returns the tid of the calling process.";

static PyObject *pypvm_mytid(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int tid = pvm_mytid();
    if (was_error(tid)) return NULL;
    resultobj = Py_BuildValue("i",tid);
    return resultobj;
}

static char pypvm_nrecv__doc[]=
"pypvm.nrecv(tid,msgtag) is a non-blocking receive.\n"
"\n"
"The  routine  pypvm.nrecv  checks  to  see if a message with\n"
"label msgtag has arrived from tid.  and  also  clears  the\n"
"current  receive  buffer if any, If a matching message has\n"
"arrived pypvm.nrecv immediately places the message in a  new\n"
"active  receive  buffer, and returns the buffer identifier\n"
"in bufid.\n"
"\n"
"If the requested message has not arrived,  then  pypvm.nrecv\n"
"immediately  returns  with  a  0  in bufid. \n"
"\n"
"A -1 or a missing in msgtag or tid matches anything. \n"
;

static PyObject *pypvm_nrecv(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  bufid;
    int  tid=-1;
    int  msgtag=-1;
    static char *kwlist[] = { "tid", "msgtag" , NULL };

    self = self;
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"|ii:pvm_nrecv", kwlist,&tid,&msgtag)) return NULL;
    bufid = pvm_nrecv(tid,msgtag);
    if (was_error(bufid)) return NULL;
    resultobj = Py_BuildValue("i",bufid);
    return resultobj;
}

static char pypvm_notify__doc[]=
"pypvm.notify(what,msgtag,[tids],[count]) requests PVM to notify the caller on\n"
"detecting certain events.  One or more notify messages (see below) are\n"
"sent by PVM back to the calling task.  The messages have tag [msgtag]\n"
"supplied to notify.\n"
"\n"
"[what] is the type  of  event  to  trigger   the   notification.\n"
"Presently one of:\n"
" pypvm.notifyDict[TaskExit]   =>      Task exits or is killed\n"
" pypvm.notifyDict[HostDelete] =>      Host is deleted or crashes\n"
" pypvm.notifyDict[HostAdd]    =>      New host is added\n"
"\n"
"[msgtag]  Message tag to be used in notification.\n"
"\n"
"[tids] For TaskExit and HostDelete, of list of task ov pvmd TIDs\n"
"to be notified about.  It is not used when [what] is HostAdd.\n"
"\n"
"[count] for HostAdd,  it determines how many messages will be sent.\n"
"It is not used for TaskExit or HostDelete.\n"
"\n"
"The notification messages have the following format:\n"
"\n"
"pypvm.notifyDict[TaskExit]\n"
"  One notify message for  each  TID  requested.   The\n"
"  message  body contains a single TID of exited task.\n"
"\n"
"pypvm.notifyDict[HostDelete]\n"
"  One notify message for  each  TID  requested.   The message  body \n"
"  contains a single pvmd-TID of exited pvmd.\n"
"\n"
"pypvm.notifyDict[HostAdd]\n"
" count notify messages are sent, one each time the local pvmds host\n"
" table is updated.  The message body contains an integer length\n"
" followed by a list of pvmd-TIDs of new pvmds.  The counter of HostAdd\n"
" messages yet to be sent is replaced by successive calls to pvm_notify.\n"
" Specifying a count of -1 turns on PvmHostAdd messages until a future\n"
" notify; a count of zero disables them.  count defaults to 0.\n"
"\n"
"TIDs in the notify messages are packed as integers.  The calling task\n"
"is responsible for receiving messages with the specified tag and\n"
"taking appropriate action.  Future versions of PVM may expand the list\n"
"of available notification events.\n"
"";
 /* 
  */

static PyObject * pypvm_notify(PyObject * self, PyObject * args, PyObject * keywords) {
  int what, msgtag, cnt=0, *tids, i;
  int info;
  PyObject *tid_list, *item;
  char * kwlist[] = {"what", "msgtag", "tids", "count", NULL };
  if (!PyArg_ParseTupleAndKeywords(args, keywords,"ii|Oi", kwlist, &what, &msgtag, &tid_list, &cnt))
     return NULL;
  if (what == PvmHostAdd) {
    /* tids not used for PvmHostAdd */
    info = pvm_notify(what,msgtag,cnt, NULL);
    return_none(info);
  }
  if (what != PvmHostAdd) {
    cnt = PyObject_Length(tid_list);
    if (cnt < 0) { 
      PyErr_SetString(PyExc_TypeError,"[tids] argument ... expected list of tids");
      return NULL;
    }
    if ((tids = (int *) PyMem_Malloc(cnt * sizeof(int))) == NULL)
       return NULL;
    for (i = 0; i < cnt; i++) {
      item = PyList_GetItem(tid_list, i);
      if (!PyInt_Check(item)) {
	/* ASSERT: second argument was not a list */
	PyErr_SetString(PyExc_TypeError,"[tids] argument: expected list of tids");
	return NULL;
      }
      tids[i] = PyInt_AsLong(item);
    }
    info = pvm_notify(what, msgtag, cnt, tids);
    PyMem_Free(tids);
    return_none(info);
    }

}

static char pypvm_parent__doc[]=
"\n"
"pypvm.parent() returns the tid of the process that spawned the calling process.    \n"
;

static PyObject *pypvm_parent(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  tid;

    tid =pvm_parent();
    if (was_error(tid)) return NULL;
    resultobj = Py_BuildValue("i",tid);
    return resultobj;
}

static char pypvm_perror__doc[]=
"\n"
"pvm_perror(msg) returns the error message of the last PVM call. The user can\n"
"(optionally) use msg to add additional information  to  the error message, \n"
"for example, its location.  I cant imagine this function being used all that\n"
"often.\n"
;

/* Hmm... should this function be used as the detail field for raised exceptions? */
static PyObject *pypvm_perror(PyObject *self, PyObject *args, PyObject * keywords) {
    char * msg = "";
    static char *kwlist[] = { "msg" , NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"|s:pvm_perror", kwlist,&msg)) {   return NULL;  }
    pvm_perror(msg);
    Py_INCREF(Py_None);
    return Py_None;
}

/* this is a utility function used by pypvm */
static PyObject * pack_one_object(PyObject *object) {
  int info;
  char * string;
  long integer;
  double real;
  if (PyString_Check(object)) {
    string = PyString_AsString(object);
    info = pvm_pkstr(string);
  } else if (PyInt_Check(object)) {
    integer = PyInt_AsLong(object);
    info = pvm_pklong(&integer,1,1);
  } else if (PyFloat_Check(object)) {
    real = PyFloat_AsDouble(object);
    info = pvm_pkdouble(&real,1,1);
  } else {
    fprintf(stderr,"Gotta make this an exception as well\n");
    exit(1);
    return NULL;
  }
  return_none(info);
}

static char pypvm_pack_by_type__doc[]=
"\n"
"pypvm.pack_by_type(arg) packs the active message buffer with the\n"
"contents of [arg].  [arg] will probably be a list, although it can be\n"
"a single value.  At the moment, only ints, floats and strings can be\n"
"packed (they get packed with pvm_pkint, pvm_pkdbl and pvm_pkstr\n"
"respectively - the programmer has no control over this... future\n"
"versions will see improvements)\n"
;
/* The Python function defined by pypvm_pack takes one argument,
   hopefully a list. */
static PyObject *pypvm_pack_by_type(PyObject *self,PyObject *args, PyObject * keywords) {
  PyObject * theList, *thisElement, *returnObject;
  int size, index;

  
  if (!PyTuple_Check(args)) { fprintf(stderr,"Seriously confused.\n"); exit(1); }
  size  = PyTuple_Size(args);
  if (size < 1) { fprintf(stderr,"Gotta make this into an exception\n"); exit(1); }
  theList = PyTuple_GetItem(args,0);
  if (!PyList_Check(theList)) { 
    /* It's not a list... perhaps they just want to pack one value.  Fair enough.
       Of course,  calling the thing they gave us "theList" is terribly confusing. */
    return (pack_one_object(theList));
  }

  /* Good.  Our first argument is a list.  Let's go through it item by item. */
  size = PyList_Size(theList);
  for (index=0;index<size;index++) {
    thisElement = PyList_GetItem(theList,index);
    returnObject = pack_one_object(thisElement);
    if (returnObject == NULL) return NULL;  /* if we couldn't pack it,  error out */
    Py_DECREF(Py_None);  /* because return_none inside pack_one_object does a Py_INCREF */
  }
  Py_INCREF(Py_None);  /* we are returning Py_None, so re-increment it back up again! */
  return Py_None;
}



static char pypvm_pk__doc[]=
"\n"
"The family of functions pypvm.pk*(list,[stride]) take a list ([list])\n"
"and sends them as an array of whatever type is given in the argument\n"
"(byte, short, int, long, float, or double).  [stride] defaults to \n"
"1,  and is the separation to use.\n"
"\n"
"One day in the future (not now) these functions will have one-argument\n"
"forms so that if the argument given is not a list,  that it will pack\n"
"it as if it were a one-element list.\n"
;
/* " */


static PyObject * pypvm_pkbyte(PyObject* self, PyObject* args, PyObject *keywords) {
  int info;
  char *chars;
  char *string_deref;
  int stride=1, list_len, i;
  PyObject *char_list, *item;
  char * kwlist[] = {"list","stride", NULL };

  if (!PyArg_ParseTupleAndKeywords(args,keywords ,"O|i", kwlist, &char_list, &stride)) {
    /* ASSERT: incorrect arguments */
    return (NULL);
  }


  if ((list_len = PyList_Size(char_list)) < 0) {
	/* ASSERT: second argument was not a list */
    PyErr_SetString(PyExc_TypeError,"argument 2: expected list of chars");
    return (NULL);
  }

  chars = (char *) PyMem_Malloc(list_len * sizeof(char));
  if (NULL == chars) {   return (NULL);   }

  for (i = 0; i < list_len; i++) {
    /* extract chars from list */
    item = PyList_GetItem(char_list, i);
    if ((NULL == item) || (PyString_Size(item) > 1)) {
      /* not a list of strings or string is more than one character */
      PyErr_SetString(PyExc_TypeError,"argument 1: expected list");
      PyMem_Free(chars);
      return NULL;
    }

    string_deref= PyString_AsString(item);
    if (NULL == string_deref) {
      PyErr_SetString(PyExc_TypeError,"argument 1: expected list only of strings");
      PyMem_Free(chars);
      return NULL;
    }
    
    /* if it is not null,  then Python promises us it must be valid and safe 
       to dereference */
    chars[i] = string_deref[0];
  }
  info = pvm_pkbyte(chars, list_len, stride);
  PyMem_Free(chars);
  return_none(info);
}


static PyObject * pypvm_pkdouble(PyObject* self, PyObject* args, PyObject *keywords) {
  int stride=1, list_len, i;
  double *doubles;
  PyObject *double_list, *item;
  int info;
  char * kwlist[] = {"list","stride",NULL };
  if (!PyArg_ParseTupleAndKeywords(args, keywords, "O|i", kwlist, &double_list, &stride))
	/* ASSERT: incorrect arguments */
	return (NULL);
  if ((list_len = PyList_Size(double_list)) < 0) {
	/* ASSERT: second argument was not a list */
    PyErr_SetString(PyExc_TypeError,"argument 1: expected list of doubles");
    return (NULL);
  }
  if ((doubles = (double *) PyMem_Malloc(list_len * sizeof(double))) ==
      NULL) return (NULL);
  for (i = 0; i < list_len; i++) {
    /* extract doubles from list */
    item = PyList_GetItem(double_list, i);
    if (PyFloat_Check(item))
      doubles[i] = PyFloat_AsDouble(item);
    else if (PyInt_Check(item))
      doubles[i] = (double) PyInt_AS_LONG(item);
    else {
      PyErr_SetString(PyExc_TypeError,
		      "argument 2: expected list of doubles");
      PyMem_Free(doubles);
      return NULL;
    }
  }
  info =  pvm_pkdouble(doubles, list_len, stride);
  PyMem_Free(doubles);
  return_none(info);
}





static PyObject * pypvm_pkfloat(PyObject* self, PyObject* args, PyObject *keywords) {
  int stride=1, list_len, i;
  float *floats;
  PyObject *float_list, *item;
  int info;
  char * kwlist[] = {"list","stride",NULL };
  if (!PyArg_ParseTupleAndKeywords(args, keywords, "O|i", kwlist, &float_list, &stride))
	/* ASSERT: incorrect arguments */
	return (NULL);
  if ((list_len = PyList_Size(float_list)) < 0) {
	/* ASSERT: second argument was not a list */
    PyErr_SetString(PyExc_TypeError,"argument 1: expected list of doubles");
    return (NULL);
  }
  if ((floats = (float *) PyMem_Malloc(list_len * sizeof(float))) ==
      NULL) return (NULL);
  for (i = 0; i < list_len; i++) {
    /* extract floats from list */
    item = PyList_GetItem(float_list, i);
    if (PyFloat_Check(item))
      floats[i] = (float) PyFloat_AsDouble(item);
    else if (PyInt_Check(item))
      floats[i] = (float) PyInt_AS_LONG(item);
    else {
      PyErr_SetString(PyExc_TypeError,"argument 2: expected list of doubles");
      PyMem_Free(floats);
      return NULL;
    }
  }
  info =  pvm_pkfloat(floats, list_len, stride);
  PyMem_Free(floats);
  return_none(info);
}

static PyObject * pypvm_pkint(PyObject *self, PyObject *args, PyObject *keywords) {
    int stride=1, list_len, *ints, i;
    int info;
    char * kwlist[] = {"list", "stride", NULL };
    PyObject *int_list, *item;
    if (!PyArg_ParseTupleAndKeywords(args, keywords, "O|i", kwlist, &int_list, &stride))
	/* ASSERT: incorrect arguments */
	return (NULL);
    if ((list_len = PyList_Size(int_list)) < 0) {
	/* ASSERT: second argument was not a list */
	PyErr_SetString(PyExc_TypeError,"argument 1: expected list of ints");
	return (NULL);
    }
    if ((ints = (int *) PyMem_Malloc(list_len * sizeof(int))) == NULL)
	return (NULL);
    for (i = 0; i < list_len; i++) {
	/* extract integers from list */
	item = PyList_GetItem(int_list, i);
	if (!PyInt_Check(item)) {
	    PyErr_SetString(PyExc_TypeError,"argument 1: expected list of ints");
	    PyMem_Free(ints);
	    return NULL;
	}
	ints[i] = (int) PyInt_AsLong(item);
    }
    info = pvm_pkint(ints, list_len, stride);
    PyMem_Free(ints);
    return_none(info);
}


static PyObject* pypvm_pklong(PyObject * self, PyObject *args, PyObject *keywords) {
  int stride=1, list_len, i;
  int info;
  long *longs;
  PyObject *long_list, *item;
  char * kwlist[] = {"list","stride", NULL };
  if (!PyArg_ParseTupleAndKeywords(args, keywords, "O|i", kwlist, &long_list, &stride))
    /* ASSERT: incorrect arguments */
    return (NULL);
  if ((list_len = PyList_Size(long_list)) < 0) {
    /* ASSERT: second argument was not a list */
    PyErr_SetString(PyExc_TypeError,"argument 1: expected list of longs");
    return (NULL);
  }
  if ((longs = (long *) PyMem_Malloc(list_len * sizeof(long))) == NULL)
    return (NULL);
  for (i = 0; i < list_len; i++) {
    /* extract longs from list */
    item = PyList_GetItem(long_list, i);
    if (!PyInt_Check(item)) {
      PyErr_SetString(PyExc_TypeError,"argument 1: expected list of longs");
      PyMem_Free(longs);
      return NULL;
    }
    longs[i] = PyInt_AsLong(item);
  }
  info = pvm_pklong(longs, list_len, stride);
  PyMem_Free(longs);  
  return_none(info);
}


static PyObject* pypvm_pkshort(PyObject *self, PyObject *args, PyObject *keywords) {
  int info;
  short *shorts;
  int stride, list_len, i;
  PyObject *short_list, *item;
  char * kwlist[]  = {"list","stride", NULL };
  if (!PyArg_ParseTupleAndKeywords(args, keywords, "O|i", kwlist, &short_list, &stride))
    /* ASSERT: incorrect arguments */
    return (NULL);
  if ((list_len = PyList_Size(short_list)) < 0) {
    /* ASSERT: second argument was not a list */
    PyErr_SetString(PyExc_TypeError,"argument 1: expected list of shorts");
    return (NULL);
  }
  if ((shorts = (short *) PyMem_Malloc(list_len * sizeof(short))) == NULL)
    return (NULL);
  for (i = 0; i < list_len; i++) {
    /* extract shorts from list */
    item = PyList_GetItem(short_list, i);
    if (!PyInt_Check(item)) {
      PyErr_SetString(PyExc_TypeError,"argument 1: expected list of shorts");
      PyMem_Free(shorts);
      return NULL;
    }
    shorts[i] = PyInt_AsLong(item);
  }
  info = pvm_pkshort(shorts, list_len, stride);
  PyMem_Free(shorts);
  return_none(info);
}


static char pypvm_pkstr__doc[] =
"pypvm.pkstr(str) packs the string [str]\n"
;

static PyObject * pypvm_pkstr(PyObject * self, PyObject * args, PyObject * keywords) {
  char *str;
  int info;
  char * kwlist[] = { "str" , NULL };
  if (!PyArg_ParseTupleAndKeywords(args, keywords, "s", kwlist, &str))
    /* ASSERT: incorrect arguments */
    return (NULL);
  info = pvm_pkstr(str);
  return_none(info);
}



static char pypvm_probe__doc[]=
"\n"
"pypvm.probe(tid,msgtag) checks  to  see if a message with label msgtag has \n"
"arrived from tid.  If a matching  message has  arrived  pypvm.probe  \n"
"returns a buffer identifier.  This bufid can be used in a  pypvm.bufinfo  \n"
"call  to determine information about the message such as its source\n"
"       and length.                             \n"
;

static PyObject *pypvm_probe(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  bufid;
    int  tid=-1;
    int  msgtag=-1;
    static char *kwlist[] = { "tid", "msgtag" , NULL };

    self = self;
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"|ii:pvm_nrecv", kwlist,&tid,&msgtag)) {
      return NULL;
    }
    bufid = pvm_probe(tid,msgtag);
    if (was_error(bufid)) return NULL;
    resultobj = Py_BuildValue("i",bufid);
    return resultobj;
}

static char pypvm_pstat__doc[]=
"\n"
"pypvm.pstat(tid) returns the status of the process identified by [tid].\n"
"Also  note that pypvm.notify() can be used to notify the caller that a task \n"
"has failed.  \n"
;

static PyObject *pypvm_pstat(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  status;
    int  tid;
    static char *kwlist[] = { "tid" , NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"i:pvm_pstat", kwlist,&tid)) 
        return NULL;
    status = pvm_pstat(tid);
    resultobj = Py_BuildValue("i",status);
    return resultobj;
}

static char pypvm_recv__doc[]=
"pypvm.recv(tid,msgtag) blocks the process until a message with label \n"
"[msgtag] has arrived from [tid].  pypvm.recv() then places the message \n"
"in a new active receive buffer, which also clears the current receive buffer.\n"
"[tid] and [msgtag] can be -1,  meaning \"match anything\", which is the default\n"
"if either is not mentioned.\n"
;

static PyObject *pypvm_recv(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject * resultobj;
    int  bufid;
    int  tid=-1;
    int  msgtag=-1;
 
    static char *kwlist[] = {"tid","msgtag",NULL};

    if(!PyArg_ParseTupleAndKeywords(args,kwds,"|ii:pvm_recv", kwlist,
                                    &tid,&msgtag))  return NULL;
    bufid = pvm_recv(tid,msgtag);
    if (was_error(bufid)) return NULL;
    resultobj = Py_BuildValue("i",bufid);
    return resultobj;
}

static char pypvm_reg_tasker__doc[]=
"Registers the calling task as a PVM task starter. This function is for \n"
"folks who are writing stuff like debugger servers and so on. \n"
;

static PyObject *pypvm_reg_tasker(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  cc;

    cc = pvm_reg_tasker();
    return_none(cc); 
}

static char pypvm_send__doc[]=
"\n"
"pypvm.send(tid,msgtag) sends a message stored in the active send buffer\n"
"to the PVM process identified by [tid].  [msgtag] is  used to label the \n"
"content of the message.\n"
;

static PyObject *pypvm_send(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  info;
    int  tid;
    int  msgtag;
    static char *kwlist[] = { "tid" , "msgtag" , NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"ii:pvm_send", kwlist,&tid,&msgtag)) return NULL;
    info= pvm_send(tid,msgtag);
    return_none(info);
}

static PyObject *pypvm_psend_str(PyObject *self, PyObject *args) {
  int info, tid, msg_tag, encoding;
  char *s;

  if (!PyArg_ParseTuple(args, "iiis", &encoding, &tid, &msg_tag, &s))
      return NULL;

  info = pvm_initsend(encoding);

  if (was_error(info)) return NULL;

  info = pvm_pkstr(s);

  if (was_error(info)) return NULL;

  info = pvm_send(tid, msg_tag);

  return_none(info);
}

static char pypvm_sendsig__doc[]=
"pypvm.sendsig(tid,signum) sends the signal number [signum] to the PVM\n"
"process identified by [tid]. \n"
;

static PyObject *pypvm_sendsig(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  info;
    int  tid;
    int  signum;
    static char *kwlist[] = { "tid" , "signum" , NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"ii:pvm_sendsig", kwlist,&tid,&signum)) return NULL;
    info = pvm_sendsig(tid,signum);
    return_none(info);
}


#ifdef PVM33COMPAT
static char pypvm_setmwid__doc[]=
"pypvm.setmwid(bufid,waitid) assigns a new wait ID to a message buffer.";

static PyObject *pypvm_setmwid(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  info;
    int  bufid;
    int  waitid;
    static char *kwlist[] = { "bufid", "waitid" , NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"ii:pvm_setmwid", kwlist,&bufid,&waitid)) 
        return NULL;
    info = pvm_setmwid(bufid,waitid);
    return_none(info);
}
#endif


static char pypvm_setopt__doc[]=
"pypvm.setopt(what,value) is a general purpose function used\n"
"to set miscellaneous options in the PVM library.  The previous\n"
"value is returned.\n"
"\n"
"[what] is an integer defining what to set.  One of:\n"
"  pypvm.opt[Route]            1    Message routing policy\n"
"  pypvm.opt[DebugMask]        2    Libpvm debug mask\n"
"  pypvm.opt[AutoErr]          3    Auto error reporting\n"
"  pypvm.opt[OutputTid]        4    Stdout destination for children\n"
"  pypvm.opt[OutputCode]       5    Output message tag\n"
"  pypvm.opt[TraceTid]         6    Trace data destination for children\n"
"  pypvm.opt[TraceCode]        7    Trace message tag\n"
"  pypvm.opt[FragSize]         8    Message fragment size\n"
"  pypvm.opt[ResvTids]         9    Allow messages to reserved tags and TIDs\n"
"  pypvm.opt[SelfOutputTid]   10    Stdout destination\n"
"  pypvm.opt[SelfOutputCode]  11    Output message tag\n"
"  pypvm.opt[SelfTraceTid]    12    Trace data destination\n"
"  pypvm.opt[SelfTraceCode]   13    Trace message tag\n"
"  pypvm.opt[ShowTids]        14    pvm_catchout prints task ids with output\n"
"  pypvm.opt[PollType]        15    Message wait policy (shared memory)\n"
"  pypvm.opt[PollTime]        16    Message spinwait duration\n"
"\n"
"[val] Integer specifying new setting of option. Usually this is just any\n"
"old number,  but there are some predefined values for pypvm.opt[Route]:\n"
"  pypvm.opt[DontRoute]        1    Dont request or grant connections\n"
"  pypvm.opt[AllowDirect]      2    (Default) Dont request but allow\n"
"  pypvm.opt[RouteDirect]      3    Request and allow connections\n"
;
/* ' " */
static PyObject *pypvm_setopt(PyObject *self, PyObject *args, PyObject *keywords) {  
  int what, val, oldval;
  char * kwlist[] = {"what", "val", NULL };
  if (!PyArg_ParseTupleAndKeywords(args, keywords, "ii", kwlist, &what, &val))
     return NULL;
  oldval = pvm_setopt(what, val);
  if (was_error(oldval)) { return NULL; }
  return (PyLong_FromLong(oldval));
}


static char pypvm_setrbuf__doc[]=
"pypvm.setrbuf(bufid) switches the active receive buffer to bufid and returns \n"
"the previous buffer id.";

static PyObject *pypvm_setrbuf(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  oldbufid;
    int  bufid;
    static char *kwlist[] = { "bufid" , NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"i:pvm_setrbuf", kwlist,&bufid)) 
        return NULL;
    oldbufid = pvm_setrbuf(bufid);
    if (was_error(oldbufid)) return NULL;
    resultobj = Py_BuildValue("i",oldbufid);
    return resultobj;
}

static char pypvm_setsbuf__doc[]=
"pypvm.setsbuf(bufid) switches the active send buffer to bufid and returns \n"
"the previous buffer id.";

static PyObject *pypvm_setsbuf(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  oldbufid;
    int  bufid;
    static char *kwlist[] = { "bufid", NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"i:pvm_setrbuf", kwlist,&bufid)) 
        return NULL;
    oldbufid = pvm_setsbuf(bufid);
    if (was_error(oldbufid)) return NULL;
    resultobj = Py_BuildValue("i",oldbufid);
    return resultobj;
}

static char pypvm_spawn__doc[]=
"pypvm.spawn(task,argv,flag, where, ntasks) starts [ntask] copies of the \n"
"executable named [task].  The tids are returned in a list.\n"
;

static PyObject *pypvm_spawn(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  numt;
    char *task;
    PyObject *argvObject;
    char **argv;
    int flag;  
    char * where;
    int ntasks;
    int * tids;
    int size,index;
    static char *kwlist[] = { "task" , "argv", "flag", "where", "ntasks" , NULL };

    if (!PyArg_ParseTupleAndKeywords(args,keywords,"sOisi", kwlist,&task,&argvObject,&flag,&where,&ntasks)) {
      /* Hmm... I'd much rather have flags=default as default */
      return NULL;
    }

    size = PyList_Size(argvObject);
    argv = calloc(size+1,sizeof(char *));
    for (index=0;index<size;index++) {
      argv[index]=PyString_AsString(PyList_GetItem(argvObject,index));
    }
    argv[size]=NULL;

    tids=calloc(ntasks,sizeof(int));
    numt = pvm_spawn(task,argv,flag,where,ntasks,tids);
    free(argv);
    if (was_error(numt)) {
      free(tids);
      return NULL;
    }
    
    resultobj = PyList_New(ntasks);
    for (index=0;index<ntasks;index++) {
      PyList_SetItem(resultobj,index,PyInt_FromLong((long) tids[index]));
    }

    free(tids);
    return resultobj;
}

/* Is this one relevant? */
/*
static PyObject *pypvm_start_pvmd(PyObject *self, PyObject *args) {
    PyObject * resultobj;
    int  result;
    int  arg0;
    char ** arg1;
    int  arg2;
    char * argc1 = 0;

    self = self;
    if(!PyArg_ParseTupleAndKeywords(args,keywords,"isi:pvm_start_pvmd", kwlist,&arg0,&argc1,&arg2)) 
        return NULL;
    if (argc1) {
        if (SWIG_GetPtr(argc1,(void **) &arg1,"_char_pp")) {
            PyErr_SetString(PyExc_TypeError,"Type error in argument 2 of pvm_start_pvmd. Expected _char_pp.");
        return NULL;
        }
    }
    result = (int )pvm_start_pvmd(arg0,arg1,arg2);
    resultobj = Py_BuildValue("i",result);
    return resultobj;
}
*/

static char pypvm_tidtohost__doc[]=
"pypvm.tidtohost(tid) returns the tid of the daemon process on the host running\n"
"the process identified by [tid] is located.\n"
;

static PyObject *pypvm_tidtohost(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    int  dtid;
    int  tid;
    static char *kwlist[] = {"tid" , NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"i:pvm_tidtohost", kwlist,&tid)) return NULL;
    dtid = pvm_tidtohost(tid);
    if (was_error(dtid)) return NULL;
    resultobj = Py_BuildValue("i",dtid);
    return resultobj;
}

char pypvm_trecv__doc[] =
"pypvm.trecv(timeout,msgtag,tid) blocks the process until a message with\n"
"label [msgtag] has arrived from [tid].  pypvm.trecv then places the\n"
"message in a new active receive buffer, also clearing the current\n"
"receive buffer.  If no matching message arrives within the specified\n"
"waiting time [timeout] (in seconds), pypvm.trecv returns without a\n"
"message.\n"
"\n"
"As usual [tid] and [msgtag] can be -1,  which is the default if either\n"
"are missing.\n"
"\n"
"pypvm.trecv returns the bufid of the new active receive buffer.\n"
;

static PyObject *pypvm_trecv(PyObject *self, PyObject *args, PyObject *keywords) {
    int tid = -1, msgtag = -1, bufid;
    double secs;
    int secs_int;
    int usecs;
    struct timeval timeout;
    char *kwlist[] = { "timeout", "msgtag" , "tid" , NULL };

    if (!PyArg_ParseTuple(args, "d|ii",&secs,  &msgtag, &tid))
	/* ASSERT: incorrect arguments */
	return (NULL);
    secs_int = (int) secs;
    timeout.tv_sec = secs_int;
    secs -= (double) secs_int;
    usecs = (double) (secs * 1000000.0);
    timeout.tv_usec = (int) usecs;
    bufid = pvm_trecv(tid, msgtag, &timeout);
    if (was_error(bufid)) { return NULL; }
    return (PyInt_FromLong(bufid));

}

static char pypvm_unexport__doc[]=
"\n"
"pypvm.unexport(name) is provided for convenience in editing the environment  \n"
"variable PVM_EXPORT, while maintaining the colon-separated list syntax it\n"
"requires.\n"
;

static PyObject *pypvm_unexport(PyObject *self, PyObject *args, PyObject * keywords) {
    PyObject * resultobj;
    char * name;
    static char *kwlist[] = { "name" , NULL };

    if(!PyArg_ParseTupleAndKeywords(args,keywords,"s:pvm_unexport", kwlist,&name)) 
        return NULL;
    pvm_unexport(name);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *pypvm_unpack(PyObject *self, PyObject *args, PyObject *keywords) {
  /* Obviously,  I have to do something with this!!!!!!!!! */
  return_none(PvmNotImpl);
}


static char pypvm_upklong__doc[]=
"Each of the pypvm.upk*(nitem,stride) routines unpacks an array of  the\n"
"given data type from the active receive buffer.  The arguments for each \n"
"of the routines are a pointer to the  array to  be  unpacked  into, \n"
"[nitem] which is the total number of items to unpack, and the (optional)\n"
"[stride] which is the stride to use when unpacking.";
 

static PyObject *pypvm_upklong(PyObject *self, PyObject *args, PyObject *keywords) {
  PyObject * resultobj,*currentobj;
  int nitem, stride = 1;
  long * data;
  int info;
  int index;
  static char *kwlist[] = { "nitem", "stride", NULL };

  if (!PyArg_ParseTupleAndKeywords(args,keywords,"i|i:pvm_upkint", kwlist,&nitem,&stride)) return NULL;
  data = (long *) calloc(nitem,sizeof(long));
  info = pvm_upklong(data,nitem,stride);
  if (was_error(info)) { free(data); return NULL; }
  resultobj = PyList_New(nitem);
  for (index=0;index<nitem;index++) {
    currentobj = PyInt_FromLong(data[index]);
    PyList_SetItem(resultobj,index,currentobj);
  }
  free(data);
  return resultobj;  
}



static char pypvm_upkint__doc[]=
"Each of the pypvm.upk*(nitem,stride) routines unpacks an array of  the\n"
"given data type from the active receive buffer.  The arguments for each \n"
"of the routines are a pointer to the  array to  be  unpacked  into, \n"
"[nitem] which is the total number of items to unpack, and the (optional)\n"
"[stride] which is the stride to use when unpacking.";
 

static PyObject *pypvm_upkint(PyObject *self, PyObject *args, PyObject *keywords) {
  PyObject * resultobj,*currentobj;
  int nitem, stride = 1;
  int * data;
  int info;
  int index;
  static char *kwlist[] = { "nitem", "stride", NULL };

  if (!PyArg_ParseTupleAndKeywords(args,keywords,"i|i:pvm_upkint", kwlist,&nitem,&stride)) return NULL;
  data = (int *) calloc(nitem,sizeof(int));
  info = pvm_upkint(data,nitem,stride);
  if (was_error(info)) { free(data); return NULL; }
  resultobj = PyList_New(nitem);
  for (index=0;index<nitem;index++) {
    currentobj = PyInt_FromLong((long) data[index]);
    PyList_SetItem(resultobj,index,currentobj);
  }
  free(data);
  return resultobj;  
}






static char pypvm_upkbyte__doc[]=
"Each of the pypvm.upk*(nitem,stride) routines unpacks an array of  the\n"
"given data type from the active receive buffer.  The arguments for each \n"
"of the routines are a pointer to the  array to  be  unpacked  into, \n"
"[nitem] which is the total number of items to unpack, and the (optional)\n"
"[stride] which is the stride to use when unpacking.";
 

static PyObject *pypvm_upkbyte(PyObject *self, PyObject *args, PyObject *keywords) {
  PyObject * resultobj,*currentobj;
  int nitem, stride = 1;
  char * data;
  int info;
  int index;
  static char *kwlist[] = { "nitem", "stride", NULL };

  if (!PyArg_ParseTupleAndKeywords(args,keywords,"i|i:pvm_upkint", kwlist,&nitem,&stride)) return NULL;
  data = (char *) calloc(nitem,sizeof(char));
  info = pvm_upkbyte(data,nitem,stride);
  if (was_error(info)) { free(data); return NULL; }
  resultobj = PyList_New(nitem);
  for (index=0;index<nitem;index++) {
    currentobj = PyInt_FromLong((long) data[index]);
    PyList_SetItem(resultobj,index,currentobj);
  }
  free(data);
  return resultobj;  
}







static char pypvm_upkshort__doc[]=
"Each of the pypvm.upk*(nitem,stride) routines unpacks an array of  the\n"
"given data type from the active receive buffer.  The arguments for each \n"
"of the routines are a pointer to the  array to  be  unpacked  into, \n"
"[nitem] which is the total number of items to unpack, and the (optional)\n"
"[stride] which is the stride to use when unpacking.";
 

static PyObject *pypvm_upkshort(PyObject *self, PyObject *args, PyObject *keywords) {
  PyObject * resultobj,*currentobj;
  int nitem, stride = 1;
  short * data;
  int info;
  int index;
  static char *kwlist[] = { "nitem", "stride", NULL };

  if (!PyArg_ParseTupleAndKeywords(args,keywords,"i|i:pvm_upkint", kwlist,&nitem,&stride)) return NULL;
  data = (short *) calloc(nitem,sizeof(short));
  info = pvm_upkshort(data,nitem,stride);
  if (was_error(info)) { free(data); return NULL; }
  resultobj = PyList_New(nitem);
  for (index=0;index<nitem;index++) {
    currentobj = PyInt_FromLong((long) data[index]);
    PyList_SetItem(resultobj,index,currentobj);
  }
  free(data);
  return resultobj;  
}



static char pypvm_upkfloat__doc[]=
"Each of the pypvm.upk*(nitem,stride) routines unpacks an array of  the\n"
"given data type from the active receive buffer.  The arguments for each \n"
"of the routines are a pointer to the  array to  be  unpacked  into, \n"
"[nitem] which is the total number of items to unpack, and the (optional)\n"
"[stride] which is the stride to use when unpacking.";

static PyObject *pypvm_upkfloat(PyObject *self, PyObject *args, PyObject *keywords) {
  PyObject * resultobj,*currentobj;
  int nitem, stride = 1;
  float * data;
  int info;
  int index;
  static char *kwlist[] = {"nitem", "stride" , NULL };

  if (!PyArg_ParseTupleAndKeywords(args,keywords,"i|i:pvm_upkint", kwlist,&nitem,&stride)) return NULL;
  data = (float *) calloc(nitem,sizeof(float));
  info = pvm_upkfloat(data,nitem,stride);
  if (was_error(info)) { free(data); return NULL; }
  resultobj = PyList_New(nitem);
  for (index=0;index<nitem;index++) {
    currentobj = PyFloat_FromDouble((double) data[index]);
    PyList_SetItem(resultobj,index,currentobj);
  }
  free(data);
  return resultobj;  
}



static char pypvm_upkdouble__doc[]=
"Each of the pypvm.upk*(nitem,stride) routines unpacks an array of  the\n"
"given data type from the active receive buffer.  The arguments for each \n"
"of the routines are a pointer to the  array to  be  unpacked  into, \n"
"[nitem] which is the total number of items to unpack, and the (optional)\n"
"[stride] which is the stride to use when unpacking.";

static PyObject *pypvm_upkdouble(PyObject *self, PyObject *args, PyObject *keywords) {
  PyObject * resultobj,*currentobj;
  int nitem, stride = 1;
  double * data;
  int info;
  int index;
  static char *kwlist[] = {"nitem", "stride" , NULL };

  if (!PyArg_ParseTupleAndKeywords(args,keywords,"i|i:pvm_upkdouble", kwlist,&nitem,&stride)) return NULL;
  data = (double *) calloc(nitem,sizeof(double));
  info = pvm_upkdouble(data,nitem,stride);
  if (was_error(info)) { free(data); return NULL; }
  resultobj = PyList_New(nitem);
  for (index=0;index<nitem;index++) {
    currentobj = PyFloat_FromDouble(data[index]);
    PyList_SetItem(resultobj,index,currentobj);
  }
  free(data);
  return resultobj;  
}

static char pypvm_upkstr__doc[]=
"pypvm.upkstr() unpacks a string from the active receive buffer";
 /* " */

static PyObject *pypvm_upkstr(PyObject *self,PyObject *args, PyObject *keywords) {
  PyObject * resultobj;
  int info; 
  char * data; 
  int bytes, msgtag, tid; 
  char * huge_buffer;   
  int actual_length;

  int ret_info;
  info = pvm_bufinfo(pvm_getrbuf(), &bytes, &msgtag, &tid);
  if (was_error(info)) { return NULL; }
  /* mallocs size of entire buffer, not necessarily just string */
  /* This is not a problem since PyString_FromString copies into a
     perfectly sized buffer (as long as the string is NULL terminated,
     which I think PVM guarantees). */
  
  huge_buffer = (char *) PyMem_Malloc(bytes * sizeof(char) + 1); 
  if (huge_buffer == NULL) { return NULL; }
  ret_info = pvm_upkstr(huge_buffer); 
  if (was_error(ret_info)) { PyMem_Free(huge_buffer); return NULL; }
  resultobj = PyString_FromString(huge_buffer);
  PyMem_Free(huge_buffer);
  return resultobj;
}


char pypvm_version__doc[] =
"\n"
"Returns the version string of the PVM this was built against.\n"
;

static PyObject *pypvm_version(PyObject *self, PyObject *args, PyObject *keywords) {
    PyObject * resultobj;
    char * result;

    result = pvm_version();
    resultobj = Py_BuildValue("s", result);
    return resultobj;
}



/* ********************************************************************** */
/*   Context manipulation functions                                       */
/* ********************************************************************** */

static char pypvm_newcontext__doc[] =
"pypvm.newcontext() returns a newly allocated context.  However, this new\n"
"context is not yet active.  See pypvm.setcontext.\n"
;
static PyObject *pypvm_newcontext(PyObject *self, PyObject *args, PyObject *keywords) {
  int info;
  info = pvm_newcontext();
  if (was_error(info)) { return NULL; }
  return Py_BuildValue("i", info);
};



static char pypvm_setcontext__doc[] =
"pypvm.setcontext(newctx) changes the current context from old_ctx (which is\n"
"return to [newctx].\n"
;
static PyObject *pypvm_setcontext(PyObject *self, PyObject *args, PyObject *keywords) {
  int info;
  int new_context;
  char * kwlist[] = {"newctx", NULL };
  if (!PyArg_ParseTupleAndKeywords(args,keywords,"i", kwlist, &new_context)) {
    return NULL;
  }
  info = pvm_setcontext(new_context);
  if (was_error(info)) { return NULL; }
  return Py_BuildValue("i", info);
};

static char pypvm_freecontext__doc[] =
"pypvm.freecontext(ctx) frees [ctx] so that it may be reused.  Contexts\n"
"are a system resource that will be exhausted if not recycled.\n"
;

static PyObject *pypvm_freecontext(PyObject *self, PyObject *args, PyObject *keywords) {
  int info;
  int context;
  char * kwlist[] = {"ctx", NULL };
  if (!PyArg_ParseTupleAndKeywords(args,keywords,"i", kwlist, &context)) {
    return NULL;
  }
  info = pvm_freecontext(context);
  return_none(info);
};

static char pypvm_getcontext__doc[] =
"pypvm.getcontext() returns the current context of the requesting task.\n"
;


static PyObject *pypvm_getcontext(PyObject *self, PyObject *args, PyObject *keywords) {
  int info;
  info = pvm_getcontext();
  if (was_error(info)) { return NULL; }
  return Py_BuildValue("i", info);
};

static PyMethodDef pvmMethods[] = {
  {"psend_str", (PyCFunction) pypvm_psend_str, METH_VARARGS},
  { "newcontext", (PyCFunction) pypvm_newcontext, METH_VARARGS| METH_KEYWORDS, pypvm_newcontext__doc },
  { "setcontext", (PyCFunction) pypvm_setcontext, METH_VARARGS| METH_KEYWORDS, pypvm_setcontext__doc },
  { "getcontext", (PyCFunction) pypvm_getcontext, METH_VARARGS| METH_KEYWORDS, pypvm_getcontext__doc },
  { "freecontext", (PyCFunction) pypvm_freecontext, METH_VARARGS| METH_KEYWORDS, pypvm_freecontext__doc },
 { "version", (PyCFunction) pypvm_version, METH_VARARGS|METH_KEYWORDS , pypvm_version__doc },
 { "unexport", (PyCFunction) pypvm_unexport, METH_VARARGS|METH_KEYWORDS, pypvm_unexport__doc },
 /*	 { "trecv", pypvm_trecv, METH_VARARGS|METH_KEYWORDS }, */
 { "tidtohost",(PyCFunction)  pypvm_tidtohost, METH_VARARGS|METH_KEYWORDS , pypvm_tidtohost__doc },
 /*{ "start_pvmd", pypvm_start_pvmd, METH_VARARGS|METH_KEYWORDS }, */
 { "spawn", (PyCFunction) pypvm_spawn, METH_VARARGS|METH_KEYWORDS , pypvm_spawn__doc }, 
 { "setsbuf", (PyCFunction) pypvm_setsbuf, METH_VARARGS|METH_KEYWORDS , pypvm_setsbuf__doc },
 { "setrbuf", (PyCFunction) pypvm_setrbuf, METH_VARARGS|METH_KEYWORDS , pypvm_setrbuf__doc },
 { "setopt", (PyCFunction) pypvm_setopt, METH_VARARGS|METH_KEYWORDS, pypvm_setopt__doc }, 
#ifdef PVM33COMPAT
 { "setmwid", (PyCFunction) pypvm_setmwid, METH_VARARGS|METH_KEYWORDS , pypvm_setmwid__doc },
#endif
 { "sendsig", (PyCFunction) pypvm_sendsig, METH_VARARGS|METH_KEYWORDS , pypvm_sendsig__doc },
 { "send", (PyCFunction) pypvm_send, METH_VARARGS|METH_KEYWORDS , pypvm_send__doc },
 /* { "register_tasker", (PyCFunction) pypvm_reg_tasker, METH_VARARGS|METH_KEYWORDS }, */
 { "recv", (PyCFunction) pypvm_recv, METH_VARARGS|METH_KEYWORDS, pypvm_recv__doc },
 { "pstat", (PyCFunction) pypvm_pstat, METH_VARARGS|METH_KEYWORDS ,pypvm_pstat__doc },
 { "probe", (PyCFunction) pypvm_probe, METH_VARARGS|METH_KEYWORDS , pypvm_probe__doc },
 { "perror", (PyCFunction) pypvm_perror, METH_VARARGS|METH_KEYWORDS , pypvm_perror__doc },
 { "parent", (PyCFunction) pypvm_parent, METH_VARARGS|METH_KEYWORDS , pypvm_parent__doc },
 { "nrecv", (PyCFunction) pypvm_nrecv, METH_VARARGS|METH_KEYWORDS , pypvm_nrecv__doc },
 { "mytid", (PyCFunction) pypvm_mytid, METH_VARARGS|METH_KEYWORDS , pypvm_mytid__doc },
 { "mstat", (PyCFunction) pypvm_mstat, METH_VARARGS|METH_KEYWORDS , pypvm_mstat__doc },
 { "mkbuf", (PyCFunction) pypvm_mkbuf, METH_VARARGS|METH_KEYWORDS , pypvm_mkbuf__doc },
 { "mcast", (PyCFunction) pypvm_mcast, METH_VARARGS|METH_KEYWORDS , pypvm_mcast__doc }, 
 { "lvgroup", (PyCFunction) pypvm_lvgroup, METH_VARARGS|METH_KEYWORDS , pypvm_lvgroup__doc },
 { "lookup", (PyCFunction) pypvm_lookup, METH_VARARGS|METH_KEYWORDS , pypvm_lookup__doc },
 { "kill", (PyCFunction) pypvm_kill, METH_VARARGS|METH_KEYWORDS , pypvm_kill__doc },
 { "joingroup",(PyCFunction) pypvm_joingroup,METH_VARARGS|METH_KEYWORDS,pypvm_joingroup__doc },
 { "insert", (PyCFunction) pypvm_insert, METH_VARARGS|METH_KEYWORDS , pypvm_insert__doc },
 { "initsend", (PyCFunction) pypvm_initsend, METH_VARARGS|METH_KEYWORDS , pypvm_initsend__doc },
 { "halt", (PyCFunction) pypvm_halt, METH_VARARGS|METH_KEYWORDS , pypvm_halt__doc },
 { "gsize", (PyCFunction) pypvm_gsize, METH_VARARGS|METH_KEYWORDS , pypvm_gsize__doc},
 { "gettid", (PyCFunction) pypvm_gettid, METH_VARARGS|METH_KEYWORDS , pypvm_gettid__doc },
 { "getsbuf", (PyCFunction) pypvm_getsbuf, METH_VARARGS|METH_KEYWORDS , pypvm_getsbuf__doc },
 { "getrbuf", (PyCFunction) pypvm_getrbuf, METH_VARARGS|METH_KEYWORDS , pypvm_getrbuf__doc },
 { "getopt",(PyCFunction)  pypvm_getopt, METH_VARARGS|METH_KEYWORDS , pypvm_getopt__doc }, 
#ifdef PVM33COMPAT
 { "getmwid", (PyCFunction) pypvm_getmwid, METH_VARARGS|METH_KEYWORDS, pypvm_getmwid__doc },
#endif
 { "getinst", (PyCFunction) pypvm_getinst, METH_VARARGS|METH_KEYWORDS , pypvm_getinst__doc },
 /* { "gather", (PyCFunction) pypvm_gather, METH_VARARGS|METH_KEYWORDS}, */
 /* { "getfds", (PyCFunction) pypvm_getfds, METH_VARARGS|METH_KEYWORDS }, */
 { "freebuf", (PyCFunction) pypvm_freebuf, METH_VARARGS|METH_KEYWORDS , pypvm_freebuf__doc },
 { "export", (PyCFunction) pypvm_export, METH_VARARGS|METH_KEYWORDS , pypvm_export__doc },
 { "exit", (PyCFunction) pypvm_exit, METH_VARARGS|METH_KEYWORDS, pypvm_exit__doc },
 { "delhosts", (PyCFunction) pypvm_delhosts, METH_VARARGS|METH_KEYWORDS , pypvm_delhosts__doc}, 
 { "addhosts", (PyCFunction) pypvm_addhosts, METH_VARARGS|METH_KEYWORDS , pypvm_addhosts__doc}, 
 { "delete", (PyCFunction) pypvm_delete, METH_VARARGS|METH_KEYWORDS, pypvm_delete__doc },
 { "tasks", (PyCFunction) pypvm_tasks, METH_VARARGS|METH_KEYWORDS , pypvm_tasks__doc },
 { "hostinfo", (PyCFunction) pypvm_hostinfo, METH_VARARGS|METH_KEYWORDS , pypvm_hostinfo__doc},
 { "bufinfo", (PyCFunction) pypvm_bufinfo, METH_VARARGS|METH_KEYWORDS , pypvm_bufinfo__doc },
  { "bcast", (PyCFunction) pypvm_bcast, METH_VARARGS|METH_KEYWORDS , pypvm_bcast__doc },
 { "barrier", (PyCFunction) pypvm_barrier, METH_VARARGS|METH_KEYWORDS , pypvm_barrier__doc },
 { "archcode",(PyCFunction)  pypvm_archcode, METH_VARARGS|METH_KEYWORDS , pypvm_archcode__doc },
 { "pack_by_type", (PyCFunction) pypvm_pack_by_type, METH_VARARGS|METH_KEYWORDS, pypvm_pack_by_type__doc },

 { "upklong", (PyCFunction) pypvm_upklong, METH_VARARGS|METH_KEYWORDS , pypvm_upklong__doc },
 { "upkint", (PyCFunction) pypvm_upkint, METH_VARARGS|METH_KEYWORDS , pypvm_upkint__doc },
 { "upkbyte", (PyCFunction) pypvm_upkbyte, METH_VARARGS|METH_KEYWORDS , pypvm_upkbyte__doc },
 { "upkshort", (PyCFunction) pypvm_upkshort, METH_VARARGS|METH_KEYWORDS , pypvm_upkshort__doc },
 { "upkdouble", (PyCFunction) pypvm_upkdouble, METH_VARARGS|METH_KEYWORDS , pypvm_upkdouble__doc },
 { "upkfloat", (PyCFunction) pypvm_upkfloat, METH_VARARGS|METH_KEYWORDS , pypvm_upkfloat__doc },
 { "upkstr", (PyCFunction) pypvm_upkstr, METH_VARARGS|METH_KEYWORDS, pypvm_upkstr__doc },

 { "pklong", (PyCFunction) pypvm_pklong, METH_VARARGS|METH_KEYWORDS , pypvm_pk__doc },
 { "pkint", (PyCFunction) pypvm_pkint, METH_VARARGS|METH_KEYWORDS , pypvm_pk__doc },
 { "pkbyte", (PyCFunction) pypvm_pkbyte, METH_VARARGS|METH_KEYWORDS , pypvm_pk__doc },
 { "pkshort", (PyCFunction) pypvm_pkshort, METH_VARARGS|METH_KEYWORDS , pypvm_pk__doc },
 { "pkdouble", (PyCFunction) pypvm_pkdouble, METH_VARARGS|METH_KEYWORDS , pypvm_pk__doc },
 { "pkfloat", (PyCFunction) pypvm_pkfloat, METH_VARARGS|METH_KEYWORDS , pypvm_pk__doc },
 { "pkstr", (PyCFunction) pypvm_pkstr, METH_VARARGS|METH_KEYWORDS, pypvm_pkstr__doc },


 { "narch", (PyCFunction) pypvm_narch, METH_VARARGS|METH_KEYWORDS, pypvm_narch__doc },
 { "notify", (PyCFunction) pypvm_notify, METH_VARARGS|METH_KEYWORDS, pypvm_notify__doc },
 { "trecv", (PyCFunction) pypvm_trecv, METH_VARARGS|METH_KEYWORDS, pypvm_trecv__doc },
 { "catchout", (PyCFunction) pypvm_catchout, METH_VARARGS|METH_KEYWORDS, pypvm_catchout__doc },
 { "config", (PyCFunction) pypvm_config, METH_VARARGS|METH_KEYWORDS, pypvm_config__doc },
 { NULL, NULL }
};


void initpypvm_core() 
{
  /* the modules aren't used yet.  One day they might be replacements for the other dictionaries */
  PyObject * exceptionDictionary, exceptionModule;
  PyObject * dataDictionary, dataModule;
  PyObject * spawnDictionary, spawnModule;
  PyObject * notifyDictionary, notifyModule;
  PyObject * resultDictionary, resultModule;
  PyObject * optDictionary, optModule;

  pvm_setopt(PvmAutoErr,0);
  /* make things quiet.  We raise an excpetion on every error anyway,  so
     it provides no extra error reporting than what we have already. */

  pypvm_module = Py_InitModule("pypvm_core", pvmMethods);
  pypvm_dictionary = PyModule_GetDict(pypvm_module);
  
  /* exceptionModule = Py_InitModule("exception",pvmExceptionMethods); */
  exceptionDictionary = PyDict_New();
  dataDictionary = PyDict_New();
  spawnDictionary = PyDict_New();
  notifyDictionary = PyDict_New();
  resultDictionary = PyDict_New();
  optDictionary = PyDict_New();
  
  /* I am wondering whether I should remove the clutter in the
     module name space by having Data and Task and ... dictionaries
     instead. */
  PyDict_SetItemString(pypvm_dictionary,"data",dataDictionary);
  PyDict_SetItemString(pypvm_dictionary,"spawnOpts",spawnDictionary);
  PyDict_SetItemString(pypvm_dictionary,"exception",exceptionDictionary);
  PyDict_SetItemString(pypvm_dictionary,"notifyDict",notifyDictionary);
  PyDict_SetItemString(pypvm_dictionary,"results",resultDictionary);
  PyDict_SetItemString(pypvm_dictionary,"opt",optDictionary);
  
  /* Data definitions for pvm_initsend and pvm_initbuf */
  PyDict_SetItemString(dataDictionary,"default",PyInt_FromLong((long) PvmDataDefault));
  PyDict_SetItemString(dataDictionary,"raw",    PyInt_FromLong((long) PvmDataRaw));
  /* The following isn't possible,  I suspect... */
  /*
  PyDict_SetItemString(dataDictionary,"inplace", PyInt_FromLong((long) PvmDataInPlace));
  PyDict_SetItemString(pypvm_dictionary,"DataFoo",     PyInt_FromLong((long) PvmDataFoo));
  */
  
  /* Flags for pvm_spawn and such like */
  PyDict_SetItemString(spawnDictionary,"TaskDefault", PyInt_FromLong((long) PvmTaskDefault));
  PyDict_SetItemString(spawnDictionary,"TaskHost",    PyInt_FromLong((long) PvmTaskHost));
  PyDict_SetItemString(spawnDictionary,"TaskArch",    PyInt_FromLong((long) PvmTaskArch));
  PyDict_SetItemString(spawnDictionary,"TaskDebug",   PyInt_FromLong((long) PvmTaskDebug));
  PyDict_SetItemString(spawnDictionary,"TaskTrace",   PyInt_FromLong((long) PvmTaskTrace));
  PyDict_SetItemString(spawnDictionary,"MppFront",    PyInt_FromLong((long) PvmMppFront));
  PyDict_SetItemString(spawnDictionary,"HostCompl",   PyInt_FromLong((long) PvmHostCompl));
  
  /* for pvm_notify */
  PyDict_SetItemString(notifyDictionary,"TaskExit", PyInt_FromLong((long) PvmTaskExit));
  PyDict_SetItemString(notifyDictionary,"HostDelete", PyInt_FromLong((long) PvmHostDelete));
  PyDict_SetItemString(notifyDictionary,"HostAdd",  PyInt_FromLong((long) PvmHostAdd));
  
  /* for pvm_setopt and pvm_getopt */
  /* I _really_ will wrap these as a better dictionary, somehow */
#define DefineOption(x)   PyDict_SetItemString(optDictionary,#x, PyInt_FromLong((long) Pvm ## x));
  DefineOption(Route);

    /* PvmRoute options */
  DefineOption(DontRoute);
  DefineOption(AllowDirect);
  DefineOption(RouteDirect);

  DefineOption(DebugMask);
  DefineOption(AutoErr);
  DefineOption(OutputTid);
  DefineOption(OutputCode);
  DefineOption(TraceTid);
  DefineOption(TraceCode);
  DefineOption(FragSize);
  DefineOption(ResvTids);
  DefineOption(SelfOutputTid);
  DefineOption(SelfOutputCode);
  DefineOption(SelfTraceTid);
  DefineOption(SelfTraceCode);
  DefineOption(ShowTids);
  DefineOption(PollType);
  DefineOption(PollConstant);
  DefineOption(PollSleep);
  DefineOption(PollTime);
  DefineOption(TaskSelf);
  DefineOption(TaskChild);


	 /* Error status values.  These will often be raised as exceptions */

  
  PyDict_SetItemString(pypvm_dictionary,"Ok", PyInt_FromLong((long) PvmOk));
  PyDict_SetItemString(exceptionDictionary,"Ok",PyInt_FromLong((long) PvmOk));
  
  
#define DefineException(x)  \
  Pvm ## x ## Exception = PyString_FromString(#x); \
  PyDict_SetItemString(exceptionDictionary, #x "Str",Pvm ## x ## Exception);  \
  Pvm ## x ## Number = PyInt_FromLong((long) Pvm ## x); \
  PyDict_SetItemString(exceptionDictionary, #x , Pvm ## x ## Number);   \
  PyDict_SetItemString(resultDictionary, #x , Pvm ## x ## Number);

  DefineException(BadParam); 
  DefineException(Mismatch);	 
  DefineException(Overflow);	 
  DefineException(NoData );	 
  DefineException(NoHost);	 
  DefineException(NoFile);	 
  DefineException(NoMem );	 
  DefineException(BadMsg);	 
  DefineException(SysErr);	 
  DefineException(NoBuf );	 
  DefineException(NoSuchBuf); 
  DefineException(NullGroup);
  DefineException(DupGroup); 
  DefineException(NoGroup);   
  DefineException(NotInGroup );
  DefineException(NoInst );	 
  DefineException(HostFail);  
  DefineException(NoParent );	 
  DefineException(NotImpl ); 
  DefineException(DSysErr);	 
  DefineException(BadVersion);	 
  DefineException(OutOfRes); 
  DefineException(DupHost);	 
  DefineException(CantStart);	 
  DefineException(Already ); 
  DefineException(NoTask);	 
  DefineException(NoEntry);	 
  DefineException(DupEntry);

  pypvmUnknownExceptionException = PyString_FromString("unknownPVMerrorCode");
  PyDict_SetItemString(exceptionDictionary,"unknownPVMerrorCodeStr",
		       pypvmUnknownExceptionException);


  /* For pvm_gather.  I don't think they'll be necessary. */
  /*
  PyDict_SetItemString(pypvm_dictionary,"PVM_STR", PyInt_FromLong((long) PVM_STR));
  PyDict_SetItemString(pypvm_dictionary,"PVM_BYTE", PyInt_FromLong((long) PVM_BYTE));
  PyDict_SetItemString(pypvm_dictionary,"PVM_SHORT", PyInt_FromLong((long) PVM_SHORT));
  PyDict_SetItemString(pypvm_dictionary,"PVM_INT", PyInt_FromLong((long) PVM_INT));
  PyDict_SetItemString(pypvm_dictionary,"PVM_FLOAT", PyInt_FromLong((long) PVM_FLOAT));
  PyDict_SetItemString(pypvm_dictionary,"PVM_CPLX", PyInt_FromLong((long) PVM_CPLX));
  PyDict_SetItemString(pypvm_dictionary,"PVM_DOUBLE", PyInt_FromLong((long)PVM_DOUBLE));
  PyDict_SetItemString(pypvm_dictionary,"PVM_DCPLX", PyInt_FromLong((long) PVM_DCPLX));
  PyDict_SetItemString(pypvm_dictionary,"PVM_LONG", PyInt_FromLong((long) PVM_LONG));
  PyDict_SetItemString(pypvm_dictionary,"PVM_USHORT", PyInt_FromLong((long)PVM_USHORT));
  PyDict_SetItemString(pypvm_dictionary,"PVM_UINT", PyInt_FromLong((long) PVM_UINT));
  PyDict_SetItemString(pypvm_dictionary,"PVM_ULONG", PyInt_FromLong((long) PVM_ULONG));
  */
}
