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

#ifndef NDEBUG
# include <algorithm>
# include <cassert>
# include <iterator>
# include <boost/assert.hpp>
#endif /* NDEBUG */

#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#ifndef __STDC_CONSTANT_MACROS
# define __STDC_CONSTANT_MACROS
#endif /* __STDC_CONSTANT_MACROS */

#include <boost/cstdint.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/properties.hpp>


#include "CL.hh"
namespace pmain = lz::main;

#include "types.hh"
namespace ptypes = lz;

#include "trie.hh"
namespace ptrie = lz::trie;


typedef boost::tokenizer <boost::escaped_list_separator <char> > Tokenizer;


int main(int argc, char *argv[])
{
  pmain::args_t args;
  pmain::CL cl;

  if (cl.parse(argc, argv, args) == EXIT_FAILURE) {
    return EXIT_FAILURE;
  }

  std::ifstream modelFile;

  modelFile.open(args.model.c_str());
  if (!modelFile.is_open()) {
    std::cerr << "Could not open file: " << args.model << std::endl;
    return EXIT_FAILURE;
  }

#ifndef NDEBUG
  std::cout << "Reading model file..." << std::endl;
#endif /* NDEBUG */

  std::string line;

  ptypes::trie::Graph g;
  ptypes::trie::Vertex root = boost::add_vertex(g);
  g[root].name = "epsilon";
  g[root].freq = 0;

  while (!modelFile.eof()) {
#ifndef NDEBUG
    std::cout << "Read line..." << line << std::endl;
#endif /* NDEBUG */

    std::getline(modelFile, line);

    if (line != "") {
      std::vector<std::string> seq;
      boost::split(seq, line, boost::is_any_of("<>"), boost::token_compress_on);

      // build the trie structure
      lz::trie::add(seq, g, root);
    }
  }
  // close the model file
  modelFile.close();

  // read the test sequences
  std::ifstream seqFile;
  seqFile.open(args.sequences.c_str());
  if (!seqFile.is_open()) {
    std::cerr << "Coult not open the sequence file: " << args.sequences << std::endl;
    return EXIT_FAILURE;
  }

#ifndef NDEBUG
  std::cout << "Reading sequence file..." << std::endl;
#endif /* NDEBUG */

  ptypes::sequences seqs;
  ptypes::sequence seq;
  boost::uint32_t curSequenceID = std::numeric_limits<boost::uint32_t>::max();

  while (!seqFile.eof()) {
#ifndef NDEBUG
    std::cout << "Read line..." << line << std::endl;
#endif /* NDEBUG */
    std::getline(seqFile, line);
    boost::trim(line);

    if (line != "") {
      std::vector<std::string> splitVec;
      boost::split(splitVec, line, boost::is_any_of(","), boost::token_compress_on);

      if (splitVec.size() != 3) {
        std::cerr << "Expect three elements, but got: " << line << std::endl;
        return EXIT_FAILURE;
      } else {
        boost::uint32_t sequenceID = boost::lexical_cast<boost::uint32_t>(splitVec[0]);

        if (sequenceID != curSequenceID) {
          seqs.push_back(seq);
          seq.clear();
          seq.push_back(splitVec[2]);
          curSequenceID = sequenceID;
        } else {
          seq.push_back(splitVec[2]);
        }
      }
    }
  }
  seqs.push_back(seq);
  seqFile.close();

  // compute log-loss
  boost::uint32_t T = 0;
  double ll = 0.0;
  for (auto ss = seqs.begin(); ss != seqs.end(); ++ss) {
    seq.clear();
    for (auto s = ss->begin(); s != ss->end(); ++s) {
      seq.push_back(*s);
#ifndef NDEBUG
      std::cout << "Compute pHat for: " << std::endl;
      std::copy(seq.begin(), seq.end(), std::ostream_iterator<std::string>(std::cout, " "));
      std::cout << std::endl;
#endif /* NDEBUG */
      double pHat = lz::trie::condProb(seq, g, root);
      ll += std::log2(pHat);
      T++;
    }
  }

  std::cout << "ll: " << ll/boost::lexical_cast<double>(T);

  boost::dynamic_properties dp;
  dp.property("id", get(&ptypes::trie::VertexProperties::id, g));
  dp.property("Name", get(&ptypes::trie::VertexProperties::name, g));
  dp.property("Frequency", get(&ptypes::trie::VertexProperties::freq, g));

  std::string outTreeFile = args.result + "/tree.gml";
  std::ofstream outTree(outTreeFile.c_str(), std::ios::out);
  boost::write_graphml(outTree, g, get(&ptypes::trie::VertexProperties::id, g), dp, false);
  outTree.close();

  std::string outGVFile = args.result + "/tree.dot";
  std::ofstream outGV(outGVFile.c_str(), std::ios::out);
  boost::write_graphviz(outGV, g, boost::make_label_writer(get(&ptypes::trie::VertexProperties::name, g)));
  outGV.close();

  return EXIT_SUCCESS;
}
