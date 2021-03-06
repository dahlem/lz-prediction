// Copyright (C) 2011, 2012, 2015 Dominik Dahlem <dominik.dahlem@gmail.com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifndef __LZ_TYPES_HH__
#define __LZ_TYPES_HH__

#ifndef __STDC_CONSTANT_MACROS
# define __STDC_CONSTANT_MACROS
#endif /* __STDC_CONSTANT_MACROS */


#include <deque>
#include <set>
#include <string>
#include <tuple>
#include <vector>


namespace lz
{

typedef std::tuple<double, boost::uint32_t> probability;
typedef std::set<std::string> set;
typedef std::vector<std::string> sequence;
typedef std::vector<sequence> sequences;
typedef std::deque<std::string> NGram;

}

#endif
