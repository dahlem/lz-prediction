// CopYRIGHT (C) 2011, 2012, 2015 Dominik Dahlem <dahlem@mit.edu>
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

/** @file CL.cc
 * Implementation of the command-line parsing of the main routine
 *
 * @author Dominik Dahlem
 */
#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifndef __STDC_CONSTANT_MACROS
# define __STDC_CONSTANT_MACROS
#endif /* __STDC_CONSTANT_MACROS */

#include <iostream>
#include <limits>

#include <boost/cstdint.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "CL.hh"


namespace lz
{
namespace main
{

CL::CL()
    : m_opt_desc(new po::options_description("Options"))
{
  // Declare the supported options.
  po::options_description opt_general("General Configuration");
  opt_general.add_options()
      (HELP.c_str(), "produce help message")
      (VERS.c_str(), "show the version")
      (RESULT.c_str(), po::value <std::string>()->default_value("./results"), "results directory.")
      (SEQUENCES.c_str(), po::value <std::string>()->default_value(""), "Sequences file")
      (MODEL.c_str(), po::value <std::string>()->default_value(""), "Model file")
      (ORDER.c_str(), po::value <boost::uint32_t>()->default_value(0), "Compute log-losses up to the given order.")
      (MAX_ORDER.c_str(), po::value <boost::uint32_t>()->default_value(std::numeric_limits<boost::uint32_t>::max()), "Limit the size of the dictionary.")
      (ADD_ALPHA.c_str(), po::value <bool>()->default_value(0), "Add alphabet to all nodes.")
      (PRINT_GRAPH.c_str(), po::value <bool>()->default_value(0), "Serialise the tree.")
      (ARGMAX_PROB.c_str(), po::value <bool>()->default_value(0), "When computing the log-loss choose the most likely transition probability.")
      (FREQ.c_str(), po::value <bool>()->default_value(0), "Frequencies given in model.")
      (SAMPLE_SIZE.c_str(), po::value <boost::uint32_t>()->default_value(0), "Sample size.")
      ;

  m_opt_desc->add(opt_general);
}


int CL::parse(int argc, char *argv[], args_t &p_args)
{
  po::variables_map vm;

  po::store(po::parse_command_line(argc, argv, (*m_opt_desc.get())), vm);
  po::notify(vm);

  if (vm.count(HELP)) {
    std::cout << (*m_opt_desc.get()) << std::endl;
    return EXIT_FAILURE;
  }

  if (vm.count(VERS)) {
    std::cout << PACKAGE_NAME << " " << PACKAGE_VERSION << std::endl;
    std::cout << argv[0] << std::endl;
    return EXIT_FAILURE;
  }

  if (vm.count(RESULT.c_str())) {
    p_args.result = vm[RESULT.c_str()].as <std::string>();
  }

  if (vm.count(SEQUENCES.c_str())) {
    p_args.sequences = vm[SEQUENCES.c_str()].as <std::string>();
  }

  if (vm.count(MODEL.c_str())) {
    p_args.model = vm[MODEL.c_str()].as <std::string>();
  }

  if (vm.count(ORDER.c_str())) {
    p_args.order = vm[ORDER.c_str()].as <boost::uint32_t>();
  }

  if (vm.count(MAX_ORDER.c_str())) {
    p_args.max_order = vm[MAX_ORDER.c_str()].as <boost::uint32_t>();
  }

  if (vm.count(ADD_ALPHA.c_str())) {
    p_args.add_alpha = vm[ADD_ALPHA.c_str()].as <bool>();
  }

  if (vm.count(PRINT_GRAPH.c_str())) {
    p_args.print_graph = vm[PRINT_GRAPH.c_str()].as <bool>();
  }

  if (vm.count(ARGMAX_PROB.c_str())) {
    p_args.argmax_prob = vm[ARGMAX_PROB.c_str()].as <bool>();
  }

  if (vm.count(FREQ.c_str())) {
    p_args.freq = vm[FREQ.c_str()].as <bool>();
  }

  if (vm.count(SAMPLE_SIZE.c_str())) {
    p_args.sample_size = vm[SAMPLE_SIZE.c_str()].as <boost::uint32_t>();
  }

  return verify(p_args);
}

int CL::verify(args_t &p_args)
{
  if (p_args.sequences == "") {
    std::cerr << "Error: Sequences file has to be specified." << std::endl;
    return EXIT_FAILURE;
  }

  if (p_args.model == "") {
    std::cerr << "Error: Model file has to be specified." << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

}
}
