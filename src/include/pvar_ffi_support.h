#ifndef __PVAR_FFI_SUPPORT_H__
#define __PVAR_FFI_SUPPORT_H__

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

#include "pgenlib_misc.h"

#ifdef __cplusplus
namespace plink2 {
#endif

struct RefcountedWptrStruct {
  uintptr_t ref_ct;
  uintptr_t p[];
};

typedef struct RefcountedWptrStruct RefcountedWptr;

RefcountedWptr* CreateRefcountedWptr(uintptr_t size);

void CondReleaseRefcountedWptr(RefcountedWptr** rwpp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PVAR_FFI_SUPPORT_H__
