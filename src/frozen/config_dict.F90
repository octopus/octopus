#include "global.h"

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME
#undef DICT_INCLUDE_PREFIX
#undef DICT_INCLUDE_HEADER
#undef DICT_INCLUDE_BODY

#define DICT_TEMPLATE_NAME config
#define DICT_TYPE_NAME json_object_t
#define DICT_TYPE_MODULE_NAME json_m

#define DICT_INCLUDE_PREFIX
#include "tdict.F90"
#undef DICT_INCLUDE_PREFIX

module config_dict_m

  use global_m
  use messages_m
  use profiling_m
  
  use json_m,  only: json_object_t
  use kinds_m, only: wp
  use strng_m, only: operator(==), string_hash=>strng_hash

  use strng_m, only:                   &
    string_t       => strng_t,         &
    string_init    => strng_init,      &
    string_get     => strng_get,       &
    string_tolower => strng_tolower,   &
    string_copy    => strng_copy,      &
    string_end     => strng_end

  implicit none

  private
    
  public ::                  &
    CONFIG_DICT_OK,          &
    CONFIG_DICT_KEY_ERROR,   &
    CONFIG_DICT_EMPTY_ERROR

  public ::             &
    config_dict_t,      &
    config_dict_len,    &
    config_dict_init,   &
    config_dict_next,   &
    config_dict_pop,    &
    config_dict_set,    &
    config_dict_get,    &
    config_dict_del,    &
    config_dict_extend, &
    config_dict_copy,   &
    config_dict_end

  public ::                 &
    config_dict_iterator_t

#define DICT_INCLUDE_HEADER
#include "tdict.F90"
#undef DICT_INCLUDE_HEADER

  integer, public, parameter :: CONFIG_DICT_NAME_LEN = 63

contains

#define DICT_INCLUDE_BODY
#include "tdict.F90"
#undef DICT_INCLUDE_BODY

end module config_dict_m

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME

!! Local Variables:
!! mode: f90
!! End:
