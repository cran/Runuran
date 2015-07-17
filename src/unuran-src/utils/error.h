/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

extern int unur_errno;
int unur_get_errno ( void );
void unur_reset_errno ( void );
const char *unur_get_strerror ( const int errnocode );
UNUR_ERROR_HANDLER *unur_set_error_handler( UNUR_ERROR_HANDLER *new_handler );
UNUR_ERROR_HANDLER *unur_set_error_handler_off( void );
