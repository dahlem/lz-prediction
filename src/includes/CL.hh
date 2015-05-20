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

/** @file CL.hh
 * Declaration of the command line parameters for disease network generation.
 *
 * @author Dominik Dahlem
 */
#ifndef __LZ_MAIN_CL_HH__
#define __LZ_MAIN_CL_HH__

#include <boost/scoped_ptr.hpp>

#include <boost/program_options/options_description.hpp>
namespace po = boost::program_options;


namespace lz
{
namespace main
{

/**
 * const variables specifying the allowed options.
 */
const std::string HELP = "help";
const std::string VERS = "version";

const std::string MODEL = "model";
const std::string SEQUENCES = "sequences";
const std::string RESULT = "result";


/** @struct
 * structure specifying the command line variables.
 */
struct args_t {
  std::string sequences;        /* sequences filename */
  std::string result;           /* results dir */
  std::string model;

  args_t(args_t const &args)
      : sequences(args.sequences), result(args.result), model(args.model)
  {}

  args_t()
      : sequences(""), result(""), model("")
  {}

  friend std::ostream& operator <<(std::ostream &p_os, const args_t &p_args)
  {
    p_os << "Parameters are:    " << std::endl << std::endl;
    p_os << "Sequence file:     " << p_args.sequences << std::endl
         << "Model file:        " << p_args.model << std::endl
         << "Results directory: " << p_args.result << std::endl;

    return p_os;
  }

};

/** @class CL
 * This class uses the boost program-options library to parse the command-line
 * parameters for the main routine of the discrete event simulator.
 */
class CL
{
 public:
  CL();

  /** @fn parse(int argc, char *argv[], args_t);
   * Parse the command-line parameters and store the relevant information
   * in a shared pointer of a structure.
   *
   * @param int number of command-line arguments
   * @param char** the command-line arguments
   * @param args_t a reference to the structure of the command-line
   *        arguments
   * @return either success or failure. In case of a failure then the help
   *        message was requested.
   */
  int parse(int, char **, args_t &);

  /** @fn verify(tDesArgsSP)
   * Verify the command-line parameters
   *
   * @param args_t a reference to the structure of the command-line
   *        arguments.
   */
  int verify(args_t &);

 private:

  /**
   * A scoped pointer to the description of the command-line arguments.
   */
  boost::scoped_ptr<po::options_description> m_opt_desc;
};

}
}

#endif
