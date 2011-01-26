/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_string * _unur_string_new ( void );
int _unur_string_append ( struct unur_string *string, const char *format, ... )
     ATTRIBUTE__FORMAT(2,3);
int _unur_string_appendtext ( struct unur_string *string, const char *text );
void _unur_string_free ( struct unur_string *string );
void _unur_string_clear ( struct unur_string *string );
