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
#include <deque>
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

  ptypes::set alphabet;
  ptypes::trie::Graph g;
  ptypes::trie::Vertex root = boost::add_vertex(g);
  g[root].name = "epsilon";
  g[root].freq = 0;

  std::string line;
  while (!modelFile.eof()) {
    std::getline(modelFile, line);

#ifndef NDEBUG
    std::cout << "Read line..." << line << std::endl;
#endif /* NDEBUG */

    if (line != "") {
      if (args.freq) {
        // expect the model input to consist of the patterns separated by <> and
        // a frequency for the leaf nodes separated by comma from the pattern
        std::vector<std::string> modelLine;
        boost::split(modelLine, line, boost::is_any_of(","), boost::token_compress_on);

        if (modelLine.size() == 2) {
          // the first part is the pattern
          ptypes::sequence seq;
          boost::split(seq, modelLine[0], boost::is_any_of("<>"), boost::token_compress_on);

          // add elements to alphabet
          std::copy(seq.begin(), seq.end(), std::inserter(alphabet, alphabet.end()));
          boost::uint32_t freq = boost::lexical_cast<boost::uint32_t>(modelLine[1]);

          // build the trie structure
          lz::trie::add(seq, g, root, args.max_order, freq);
        } else {
          std::cerr << "Error: Expect comma separated model line in the form a<>b,10" << std::endl;
        }
      } else {
        // expect the model input to consist of just the patterns separated by <>
        ptypes::sequence seq;
        boost::split(seq, line, boost::is_any_of("<>"), boost::token_compress_on);

        // add elements to alphabet
        std::copy(seq.begin(), seq.end(), std::inserter(alphabet, alphabet.end()));

        // build the trie structure
        lz::trie::add(seq, g, root, args.max_order, 1);
      }
    }
  }
  // close the model file
  modelFile.close();

  // add the alphabet to each node if it does not already exist
  if (args.add_alpha) {
    BOOST_FOREACH(ptypes::trie::Vertex v, (boost::vertices(g))) {
      ptypes::set children = alphabet;
      BOOST_FOREACH(ptypes::trie::Edge e, (boost::out_edges(v, g))) {
        std::string name = g[boost::target(e, g)].name;
        if (children.find(name) != children.end()) {
          children.erase(name);
        }
      }
      BOOST_FOREACH(std::string child, children) {
        boost::uint32_t numVertices = boost::num_vertices(g);
        ptypes::trie::Vertex newV = boost::add_vertex(g);
        std::pair<ptypes::trie::Edge, bool> e = boost::add_edge(v, newV, g);

        g[newV].id = numVertices + 1;
        g[newV].freq = 1;
        g[newV].name = child;
      }
    }
  } else {
    ptypes::set children = alphabet;
    BOOST_FOREACH(ptypes::trie::Edge e, (boost::out_edges(root, g))) {
      std::string name = g[boost::target(e, g)].name;
      if (children.find(name) != children.end()) {
        children.erase(name);
      }
    }
    BOOST_FOREACH(std::string child, children) {
      boost::uint32_t numVertices = boost::num_vertices(g);
      ptypes::trie::Vertex newV = boost::add_vertex(g);
      std::pair<ptypes::trie::Edge, bool> e = boost::add_edge(root, newV, g);

      g[newV].id = numVertices + 1;
      g[newV].freq = 1;
      g[newV].name = child;
    }
  }

  ptypes::trie::FreqPropertyMap vertex_freq_map = boost::get(&ptypes::trie::VertexProperties::freq, g);
  boost::depth_first_search(g, boost::visitor(ptypes::trie::dfs_frequencies<ptypes::trie::FreqPropertyMap>(vertex_freq_map)));

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
    std::getline(seqFile, line);
#ifndef NDEBUG
    std::cout << "Read line..." << line << std::endl;
#endif /* NDEBUG */

    if (line != "") {
      ptypes::sequence splitVec;
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
  std::vector<double> ll;
  std::vector<boost::uint32_t> T;
  ptypes::NGram ngram;

  ll.resize(1 + args.order);
  T.resize(1 + args.order);

  for (auto ss = seqs.begin(); ss != seqs.end(); ++ss) {
    ngram.clear();
    for (auto s = ss->begin(); s != ss->end(); ++s) {
      ngram.push_back(*s);
#ifndef NDEBUG
      std::cout << "Compute pHat for: " << std::endl;
      std::copy(ngram.begin(), ngram.end(), std::ostream_iterator<std::string>(std::cout, " "));
      std::cout << std::endl;
#endif /* NDEBUG */

      double pHat = 0.0;
      if (args.argmax_prob) {
        pHat = lz::trie::argmaxProb(ngram, g, root);
      } else {
        pHat = lz::trie::condProb(ngram, g, root);
      }
      ll[0] += std::log2(pHat);
    }
    T[0] += ss->size();
  }

  if (args.order > 0) {
    for (auto o = 1; o <= args.order; ++o) {
#ifndef NDEBUG
      std::cout << "Compute log-loss for order cutoff " << o << std::endl;
#endif /* NDEBUG */
      for (auto ss = seqs.begin(); ss != seqs.end(); ++ss) {
        ngram.clear();
        for (auto s = ss->begin(); s != ss->end(); ++s) {
          ngram.push_back(*s);
          if (ngram.size() > o) {
            ngram.pop_front();
          }
#ifndef NDEBUG
          std::cout << "Compute pHat for: " << std::endl;
          std::copy(ngram.begin(), ngram.end(), std::ostream_iterator<std::string>(std::cout, " "));
          std::cout << std::endl;
#endif /* NDEBUG */

          double pHat = 0.0;
          if (args.argmax_prob) {
            pHat = lz::trie::argmaxProb(ngram, g, root);
          } else {
            pHat = lz::trie::condProb(ngram, g, root);
          }
          ll[o] += std::log2(pHat);
        }
        T[o] += ss->size();
      }
    }
  }

  std::string outLogLossFile = args.result + "/log-loss.csv";
  std::ofstream outLogLoss(outLogLossFile.c_str(), std::ios::out);
  outLogLoss << "order,logloss,T" << std::endl;
  for (auto i = 0; i < ll.size(); ++i) {
#ifndef NDEBUG
    std::cout << i << "," << -ll[i]/boost::lexical_cast<double>(T[i]) << "," << T[i] << std::endl;
#endif /* NDEBUG */
    outLogLoss << i << "," << -ll[i]/boost::lexical_cast<double>(T[i]) << "," << T[i] << std::endl;
  }
  outLogLoss.close();

  if (args.print_graph) {
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
  }

  return EXIT_SUCCESS;
}
