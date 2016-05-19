#include "global.h"
#include "defs.h"




int myexit(int call_code)
{
  int cleanfinish;

  fprintf(stderr, "proc: %02d : Exiting cc: %d\n", myid, call_code);
  fclose(fail_file);
  fclose(log_file);
  if (myid <= 0) {
    if (1) {
      if (logfull_file)
	fclose(logfull_file);
    }
  }
  if (DOLOGSTEP) {
    if (1) {
      if (logstep_file)
	fclose(logstep_file);
    }
  }
  if (DOLOGPERF) {
    if (1) {
      if (logperf_file)
	fclose(logperf_file);
    }
  }
  if (DODTDIAG) {
    if (1) {
      if (logdt_file)
	fclose(logdt_file);
    }
  }
  if (DOFLOORDIAG >= 1) {
    if (1) {
      if (logfl_file)
	fclose(logfl_file);
    }
  }
  if (DOSPDIAG >= 1) {
    if (1) {
      if (logsp_file)
	fclose(logsp_file);
    }
  }

  if (call_code > 0) {
    fprintf(stderr,
	    "proc: %02d : Failure.  Please check failure file: cc: %d\n",
	    myid, call_code);

    if (call_code == 5) {
      cleanfinish = 1;
    } else {
      cleanfinish = 0;
    }
  } else
    cleanfinish = 1;

  if (cleanfinish) {
    fprintf(stderr,
	    "Ending Computation on proc: %02d, holding for other cpus\n",
	    myid);


    if (myid <= 0)
      fprintf(stderr, "Ended Computation on all processors\n");
  }

  exit(0);
  return (0);
}


void itoa(int x, char *p)
{
  int temp1;
  int i;
  int digits = 0;

  temp1 = x;
  if (temp1 == 0)
    digits = 1;
  else {
    for (i = 0; i < DIGILEN; i++) {
      if (temp1 == 0) {
	digits = i;
	break;
      } else
	temp1 = temp1 / 10;
    }
  }
  if (digits == 0) {
    fprintf(fail_file, "problem with itoa function: x: %d p: %s\n", x,
	    p);
    myexit(1);
  }

  for (i = 0; i < digits; i++) {
    temp1 = x / 10;
    temp1 = x - temp1 * 10;
    p[digits - i - 1] = ZER + temp1;
    x = x / 10;
  }
  p[digits] = '\0';
}


int mysys(char *com1, char *com2)
{

  pid_t pid = -1;
  int trycnt = 0;

  while (pid < 0) {
    trycnt++;
    pid = fork2();		/* forking new process */
    if (pid == 0) {		/* child do the command */
      execlp(com1, com1, com2, NULL);
      _exit(0);
    } else if (pid < 0) {	/* error creating new process */
      fprintf(fail_file, "can't create new process for: %s %s\n",
	      com1, com2);
      fprintf(fail_file, "try#: %d\n", trycnt);
      fprintf(fail_file, "errno: %d\n", errno);
    }
  }
  return (pid);
}

int fork2(void)
{
  pid_t pid;
  int status;

  if (!(pid = fork())) {
    switch (fork()) {
    case 0:
      return 0;
    case -1:
      _exit(errno);		/* assumes all errnos are <256 */
    default:
      _exit(0);
    }
  }

  if (pid < 0 || waitpid(pid, &status, 0) < 0)
    return -1;

  if (WIFEXITED(status))
    if (WEXITSTATUS(status) == 0)
      return 1;
    else
      errno = WEXITSTATUS(status);
  else
    errno = EINTR;		/* well, sort of :-) */

  return -1;
}
