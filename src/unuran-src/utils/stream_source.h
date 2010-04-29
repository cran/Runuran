/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

void _unur_log_printf (const char *genid, const char *filename, int line, const char *format, ...)
  ATTRIBUTE__FORMAT(4,5);
void _unur_log_debug (const char *format, ...)
  ATTRIBUTE__FORMAT(1,2);
int _unur_read_data (const char *file, int no_of_entries, double **array);
