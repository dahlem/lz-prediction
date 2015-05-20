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

/** @file trie.hh
 * Declaration of the trie structure to efficiently divide pattern sequences and their partial results.
 *
 * @author Dominik Dahlem
 */
#ifndef __LZ_TRIE_HH__
#define __LZ_TRIE_HH__

#ifndef __STDC_CONSTANT_MACROS
# define __STDC_CONSTANT_MACROS
#endif /* __STDC_CONSTANT_MACROS */

#include <cstring>

#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>

#include "types.hh"
namespace ptypes = lz;


namespace lz
{
namespace trie
{

struct VertexProperties
{
  boost::uint32_t id;
  std::string name;
  boost::uint32_t freq;
};


/** @typedef Graph
 * Specifies the graph to be used as a trie.
 */
typedef boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::bidirectionalS,
  boost::property <boost::vertex_index_t, int, VertexProperties> > Graph;

/** @typedef Vertex
 * Specifies the vertex descriptor of a graph
 */
typedef boost::graph_traits <Graph>::vertex_descriptor Vertex;

/** @typedef Edge
 * Specifies the edge descriptor of a graph
 */
typedef boost::graph_traits <Graph>::edge_descriptor Edge;

/** @typedef VertexIterator
 * Specifies the iterator for the vertices
 */
typedef boost::graph_traits <Graph>::vertex_iterator VertexIterator;

/** @typedef EdgeIterator
 * Specifies the iterator for the edges
 */
typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;

/** @typedef OutEdgeIterator
 * Specifies the iterator for the out degree edges
 */
typedef boost::graph_traits <Graph>::out_edge_iterator OutEdgeIterator;

/** @typedef InEdgeIterator
 * Specifies the iterator for the in degree edges
 */
typedef boost::graph_traits <Graph>::in_edge_iterator InEdgeIterator;


Vertex find(const std::string &p_symbol, Graph &p_graph, Vertex &p_vertex)
{
  BOOST_FOREACH(Edge e, (boost::out_edges(p_vertex, p_graph))) {
    bool found = false;

    if (p_graph[target(e, p_graph)].name == p_symbol) {
      found = true;
    }

    if (found) {
      return target(e, p_graph);
    }
  }

  return NULL;
}

void add(const ptypes::sequence &p_seq, Graph &p_graph, Vertex &p_root)
{
  Vertex curV = p_root;
  Vertex nextV = NULL;
  boost::uint32_t i = 0;
  boost::uint32_t numVertices = boost::num_vertices(p_graph);

  for (; i < p_seq.size(); ++i) {
    nextV = find(p_seq[i], p_graph, curV);

    // break, because we have to add nodes now
    if (nextV == NULL) {
      break;
    } else {
      p_graph[curV].freq++;
      curV = nextV;
    }
  }

  // add the nodes
  for (; i < p_seq.size(); ++i) {
    Vertex newV = boost::add_vertex(p_graph);
    std::pair<Edge, bool> e = boost::add_edge(curV, newV, p_graph);

    p_graph[newV].id = numVertices;
    p_graph[newV].freq = 1;
    p_graph[newV].name = p_seq[i];

    curV = newV;

    numVertices += 1;
  }
}

double condProb(ptypes::sequence &p_seq, Graph &p_graph, Vertex &p_root)
{
  Vertex curV = p_root, nextV;

  for (auto i = 0; i < p_seq.size(); ++i) {
#ifndef NDEBUG
    std::cout << "Walk: " << p_seq[i] << std::endl;
#endif /* NDEBUG */

    nextV = find(p_seq[i], p_graph, curV);

    if (nextV == NULL) {
#ifndef NDEBUG
      std::cout << "Restart" << std::endl;
#endif /* NDEBUG */
      curV = p_root;
      nextV = find(p_seq[i], p_graph, curV);

      if (nextV == NULL) {
        std::cerr << "ERROR: Cannot find first order!" << std::endl;
      }
    }
    curV = nextV;
  }

  // we know in a tree there is only one parent
  InEdgeIterator parent, in_end;
  boost::tie(parent, in_end) = boost::in_edges(curV, p_graph);

  double childFreq = boost::lexical_cast<double>(p_graph[curV].freq);
  double parentFreq = boost::lexical_cast<double>(p_graph[boost::source(*parent, p_graph)].freq);

#ifndef NDEBUG
  std::cout << "pHat = " << childFreq << "/" << parentFreq << " = " << childFreq/parentFreq << std::endl;
#endif /* NDEBUG */

  return childFreq/parentFreq;
}

}
}

#endif
