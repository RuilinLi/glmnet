// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.

#include "pvar_ffi_support.h"

namespace plink2 {

RefcountedWptr* CreateRefcountedWptr(uintptr_t size) {
  RefcountedWptr* rwp = S_CAST(RefcountedWptr*, malloc(sizeof(RefcountedWptr) + size * sizeof(intptr_t)));
  if (!rwp) {
    return nullptr;
  }
  rwp->ref_ct = 1;
  return rwp;
}

void CondReleaseRefcountedWptr(RefcountedWptr** rwpp) {
  RefcountedWptr* rwp = *rwpp;
  if (!rwp) {
    return;
  }
  --rwp->ref_ct;
  if (!rwp->ref_ct) {
    free(rwp);
  }
  *rwpp = nullptr;
}


}  // namespace plink2
