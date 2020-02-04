/* Copyright (c) 2000-2020 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <time.h>
#include <stdarg.h>
#include <ctype.h>
#ifdef R_UNURAN
#include <R_ext/Error.h>
#endif
#ifndef UNUR_LOG_FILE
#  define UNUR_LOG_FILE "unuran.log"
#endif
static FILE *_unur_logfile_open( void );  
static FILE *unur_stream = NULL;
static const char GENID_UNKNOWN[] = "UNURAN";
void
_unur_log_printf( const char *genid ATTRIBUTE__UNUSED,
		  const char *filename ATTRIBUTE__UNUSED, 
		  int line ATTRIBUTE__UNUSED, 
		  const char *format ATTRIBUTE__UNUSED, ... )
{
#ifdef UNUR_ENABLE_LOGGING
  va_list ap;
  if (!genid) genid = GENID_UNKNOWN;
  va_start(ap, format);
  if (!unur_stream) unur_get_stream();
  fprintf(unur_stream,"%s: %s:%d - ",genid,filename,line);
  vfprintf(unur_stream,format,ap);
  fprintf(unur_stream,"\n");
  fflush(unur_stream);   
  va_end(ap);
#else
  return;
#endif
} 
void _unur_log_debug( const char *format ATTRIBUTE__UNUSED, ... )
{
#ifdef UNUR_ENABLE_LOGGING
  va_list ap;
  va_start(ap, format);
  if (!unur_stream) unur_get_stream();
  vfprintf(unur_stream,format,ap);
  fflush(unur_stream);   
  va_end(ap);
#else
  return;
#endif
} 
FILE * 
unur_set_stream( FILE *new_stream )
{
  FILE * previous_stream;
  _unur_check_NULL( GENID_UNKNOWN,new_stream,NULL );
  previous_stream = unur_stream;
  unur_stream = new_stream;
  return previous_stream;
} 
FILE * 
unur_get_stream( void )
{
  if (unur_stream == NULL) {
    unur_stream = _unur_logfile_open();
  }
  return unur_stream;
} 
FILE *
_unur_logfile_open( void )
{
  static FILE* LOG = NULL;
  if (LOG) return LOG;  
#if defined(UNUR_ENABLE_LOGGING)
  LOG = fopen(UNUR_LOG_FILE,"w");
  if (!LOG) {
    fprintf(stderr,"Warning: cannot open logfile %s !\n",UNUR_LOG_FILE);
    fprintf(stderr,"Use STDERR instead !\n");
    LOG = stderr;
  }
  fprintf(LOG,"\nUNU.RAN - Universal Non-Uniform RANdom number generator\n\n");
  fprintf(LOG,"Version: %s\n",PACKAGE_VERSION);
  {
    time_t started;   
    if (time( &started ) != -1)
      fprintf(LOG,"%s",ctime(&started));
  }
  fprintf(LOG,"\n=======================================================\n\n");
#elif defined(R_UNURAN)
  LOG = fopen(UNUR_LOG_FILE,"w");
  if (!LOG) { error("Cannot open LOG file."); }
#else
  LOG = stderr;
#endif
  return LOG;
} 
int
_unur_read_data( const char *filename, int no_of_entries, double **ar )
{
#define LINELENGTH  1024      
  const int datasize = 1000; 
  int i, j;
  char *c;
  int memfactor = 1;
  char line[LINELENGTH];
  char *toline;
  char *chktoline;
  double *data;              
  int n_data;                
  FILE *fp;
  *ar = NULL;
  n_data = 0;
  if (datasize < no_of_entries) {
    _unur_error("read_data",UNUR_ERR_GEN_DATA,"No of entries > max datasize");
    return 0;   
  }
  data = _unur_xmalloc(memfactor * datasize * sizeof(double));
  fp = fopen(filename, "r");
  if (fp == NULL) {
    _unur_error("read_data",UNUR_ERR_GENERIC,"cannot open file");
    free(data);
    return 0; 
  }
  for ( c = fgets(line, LINELENGTH, fp), i=0;
        !feof(fp) && c;
        c = fgets(line, LINELENGTH, fp) ) {
    if (i > memfactor*datasize - no_of_entries-2){
      memfactor++;
      data = _unur_xrealloc(data, memfactor*datasize*sizeof(double));
    }
    if ( ! (isdigit(line[0]) || line[0] == '.' || line[0] == '+' 
           || line[0] == '-' ) )
      continue;
    ++n_data;
    toline = line;      
    for (j=0 ; j<no_of_entries; i++, j++){
      chktoline = toline;
      data[i] = strtod(toline, &toline); 
      if (chktoline == toline) {
	_unur_error("read_data",UNUR_ERR_GEN_DATA,"data file not valid");
	free(data);
	fclose(fp);
	return 0;    
      }  
    } 
  }    
  fclose(fp);
  data = _unur_xrealloc( data, (i+1) * sizeof(double) );
  *ar = data;
  return n_data;
#undef LINELENGTH
} 
