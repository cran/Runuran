/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef DEBUG_H_SEEN
#define DEBUG_H_SEEN
#define UNUR_DEBUG_OFF     (0u)           
#define UNUR_DEBUG_ALL     (~0u)      
#define UNUR_DEBUG_INIT    0x00000001u    
#define UNUR_DEBUG_SETUP   0x00000fffu    
#define UNUR_DEBUG_ADAPT   0x00fff000u    
#define UNUR_DEBUG_SAMPLE  0xff000000u    
int unur_set_debug( UNUR_PAR *parameters, unsigned debug );
int unur_chg_debug( UNUR_GEN *generator, unsigned debug );
int unur_set_default_debug( unsigned debug );
#endif 
