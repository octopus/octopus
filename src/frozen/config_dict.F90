#include "global.h"

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME
#undef DICT_INCLUDE_PREFIX
#undef DICT_INCLUDE_HEADER
#undef DICT_INCLUDE_BODY

#define DICT_TEMPLATE_NAME config
#define DICT_TYPE_NAME json_object_t
#define DICT_TYPE_MODULE_NAME json_oct_m
#define DICT_TYPE_EXTERNAL

module config_dict_oct_m

  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  
#define DICT_INCLUDE_PREFIX
#include "tdict_inc.F90"
#undef DICT_INCLUDE_PREFIX

  implicit none

  private
    
  public ::               &
    CONFIG_DICT_NAME_LEN

  public ::                  &
    CONFIG_DICT_OK,          &
    CONFIG_DICT_KEY_ERROR,   &
    CONFIG_DICT_EMPTY_ERROR

  public ::        &
    config_dict_t

  public ::             &
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
#include "tdict_inc.F90"
#undef DICT_INCLUDE_HEADER

contains

#define DICT_INCLUDE_BODY
#include "tdict_inc.F90"
#undef DICT_INCLUDE_BODY

end module config_dict_oct_m

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME

!! Local Variables:
!! mode: f90
!! End:
