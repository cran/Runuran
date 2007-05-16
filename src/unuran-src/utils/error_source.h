/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_ERROR_SOURCE_H_SEEN
#define UNUR_ERROR_SOURCE_H_SEEN
void _unur_error_x( const char *objid, const char *file, int line, 
		    const char *errortype, int errorcode, const char *reason );
#ifdef UNUR_COOKIES
void _unur_error_cookies( const char *file, int line, unsigned observed, unsigned expected );
#endif
#define _unur_error(genid,errorcode,reason) \
   do { \
      _unur_error_x((genid),__FILE__,__LINE__,"error",(errorcode),(reason)); \
   } while (0)
#define _unur_warning(genid,errorcode,reason) \
   do { \
      _unur_error_x((genid),__FILE__,__LINE__,"warning",(errorcode),(reason)); \
   } while (0)
#endif  
