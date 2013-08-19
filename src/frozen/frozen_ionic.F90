#include "global.h"

#define TEMPLATE_NAME frozen
#define SUBTEMPLATE_NAME fio
#include "tionic_term.F90"
#undef SUBTEMPLATE_NAME
#undef TEMPLATE_NAME

module frozen_ionic_m

  use frozen_ionic_term_m, only:                                       &
    frozen_ionic_t               => frozen_ionic_term_t,               &
    frozen_ionic_init            => frozen_ionic_term_init,            &
    frozen_ionic_update          => frozen_ionic_term_update,          &
    frozen_ionic_get             => frozen_ionic_term_get,             &
    frozen_ionic_get_energy      => frozen_ionic_term_get_energy,      &
    frozen_ionic_get_interaction => frozen_ionic_term_get_interaction, &
    frozen_ionic_copy            => frozen_ionic_term_copy,            &
    frozen_ionic_end             => frozen_ionic_term_end

   implicit none

  private
  public ::                       &
    frozen_ionic_t,               &
    frozen_ionic_init,            &
    frozen_ionic_update,          &
    frozen_ionic_get,             &
    frozen_ionic_get_energy,      &
    frozen_ionic_get_interaction, &
    frozen_ionic_copy,            &
    frozen_ionic_end
  
end module frozen_ionic_m

!! Local Variables:
!! mode: f90
!! End:
