/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef URNG_RNGSTREAMS_H_SEEN
#define URNG_RNGSTREAMS_H_SEEN
#include <RngStream.h>
UNUR_URNG *unur_urng_rngstream_new( const char *urngstr );
UNUR_URNG *unur_urng_rngstreamptr_new( RngStream rngstream );
#endif  
