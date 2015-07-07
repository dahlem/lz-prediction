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

#include <algorithm>
#include <cmath>
#include <cstring>
#include <deque>
#include <functional>
#include <numeric>
#include <tuple>

#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/numeric/conversion/cast.hpp>
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

/** @typedef FreqPropertyMap
 * Specifies the frequency property map
 */
typedef boost::property_map<Graph, boost::uint32_t VertexProperties::*>::type FreqPropertyMap;


Vertex find(const ptypes::NGram::value_type &p_symbol, Graph &p_graph, Vertex &p_vertex)
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

void add(const ptypes::sequence &p_seq, Graph &p_graph, Vertex &p_root, boost::uint32_t p_maxOrder, boost::uint32_t p_freq)
{
  Vertex curV = p_root;
  Vertex nextV = NULL;
  boost::uint32_t i = 0;
  boost::uint32_t numVertices = boost::num_vertices(p_graph);

  boost::uint32_t order = std::min(boost::numeric_cast<boost::uint32_t>(p_seq.size()), p_maxOrder);
  for (; i < order; ++i) {
    nextV = find(p_seq[i], p_graph, curV);
    // break, because we have to add nodes now
    if (nextV == NULL) {
#ifndef NDEBUG
      std::cout << "No edge to " << p_seq[i] << std::endl;
#endif /* NDEBUG */
      break;
    } else {
#ifndef NDEBUG
      std::cout << "Walk along " << p_seq[i] << std::endl;
#endif /* NDEBUG */
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

#ifndef NDEBUG
    std::cout << "Add edge between " << p_graph[curV].name << " and " << p_graph[newV].name << std::endl;
#endif /* NDEBUG */

    curV = newV;

    numVertices += 1;
  }

  // update the frequency of the last node added
  p_graph[curV].freq = p_freq;
}

bool tupleCompare(ptypes::probability left, ptypes::probability right)
{
  return std::get<0>(left)<std::get<0>(right);
}

ptypes::probability condProb(ptypes::NGram &p_seq, Graph &p_graph, Vertex &p_root)
{
  Vertex curV = p_root, nextV, prevV;
  auto i = 0;

  for (; i < p_seq.size(); ++i) {
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
        std::cerr << "ERROR: Cannot find first order for " << p_seq[i] << std::endl;
      }
    }
    prevV = curV;
    curV = nextV;
  }

  double childFreq = boost::lexical_cast<double>(p_graph[curV].freq);
  double parentFreq = boost::lexical_cast<double>(p_graph[prevV].freq);

#ifndef NDEBUG
  std::cout << "pHat = " << childFreq << "/" << parentFreq << " = " << childFreq/parentFreq << std::endl;
#endif /* NDEBUG */

  double result = childFreq/parentFreq;

  return std::make_tuple(result, i+1);
}

ptypes::probability argmaxProb(ptypes::NGram &p_seq, Graph &p_graph, Vertex &p_root)
{
  ptypes::NGram ngram = p_seq;
  std::vector<probability> probs;
  std::vector<ptypes::probability>::iterator argmax;

  auto p = condProb(ngram, p_graph, p_root);
  probs.push_back(p);

  for (auto i = ngram.size(); i > 1; --i) {
    ngram.pop_front();
    p = condProb(ngram, p_graph, p_root);
    probs.push_back(p);
  }

  argmax = std::max_element(probs.begin(), probs.end(), tupleCompare);
  return *argmax;
}



template <class VertexFreqMap>
class dfs_frequencies : public boost::default_dfs_visitor
{
 public:
  dfs_frequencies(VertexFreqMap &p_freqMap)
      : m_freqMap(p_freqMap) {
#ifndef NDEBUG
    std::cout << "dfs_frequencies" << std::endl;
#endif /* NDEBUG */
  }

  ~dfs_frequencies() {
#ifndef NDEBUG
    std::cout << "~dfs_frequencies" << std::endl;
#endif /* NDEBUG */
  }

  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex s, Graph &g) {
#ifndef NDEBUG
    std::cout << "Initialise: " << g[s].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Vertex, typename Graph>
  void start_vertex(Vertex s, Graph &g) {
#ifndef NDEBUG
    std::cout << "Start: " << g[s].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex u, Graph &g) {
#ifndef NDEBUG
    std::cout << "Discover: " << g[u].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Edge, typename Graph>
  void examine_edge(Edge e, Graph &g) {
#ifndef NDEBUG
    std::cout << "Examine edge: " << g[boost::source(e, g)].name << "-" << g[boost::target(e, g)].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Edge, typename Graph>
  void tree_edge(Edge e, Graph &g) {
#ifndef NDEBUG
    std::cout << "Tree edge: " << g[boost::source(e, g)].name << "-" << g[boost::target(e, g)].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Edge, typename Graph>
  void back_edge(Edge e, Graph &g) {
#ifndef NDEBUG
    std::cout << "Back edge: " << g[boost::source(e, g)].name << "-" << g[boost::target(e, g)].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Edge, typename Graph>
  void forward_or_cross_edge(Edge e, Graph &g) {
#ifndef NDEBUG
    std::cout << "Forward or Cross edge: " << g[boost::source(e, g)].name << "-" << g[boost::target(e, g)].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Edge, typename Graph>
  void finish_edge(Edge e, Graph &g) {
#ifndef NDEBUG
    std::cout << "Finish edge: " << g[boost::source(e, g)].name << "-" << g[boost::target(e, g)].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Vertex, typename Graph>
  void finish_vertex(Vertex u, Graph &g) {
#ifndef NDEBUG
    std::cout << "Finish vertex: " << g[u].name << std::endl;
#endif /* NDEBUG */

    BOOST_FOREACH(Edge e, boost::out_edges(u, g)) {
      m_freqMap[u] += m_freqMap[boost::target(e, g)];
    }
  }

 private:
  VertexFreqMap &m_freqMap;
};


template <class VertexFreqMap>
class dfs_entropies : public boost::default_dfs_visitor
{
 public:
  dfs_entropies(VertexFreqMap &p_freqMap, double &p_h, double &p_C, double &p_bic)
      : m_freqMap(p_freqMap), m_h(p_h), m_C(p_C), m_bic(p_bic) {
#ifndef NDEBUG
    std::cout << "dfs_entropies" << std::endl;
#endif /* NDEBUG */
  }

  ~dfs_entropies() {
#ifndef NDEBUG
    std::cout << "~dfs_entropies" << std::endl;
#endif /* NDEBUG */
  }

  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex s, Graph &g) {
#ifndef NDEBUG
    std::cout << "Initialise: " << g[s].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Vertex, typename Graph>
  void start_vertex(Vertex s, Graph &g) {
#ifndef NDEBUG
    std::cout << "Start: " << g[s].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex u, Graph &g) {
#ifndef NDEBUG
    std::cout << "Discover: " << g[u].name << std::endl;
#endif /* NDEBUG */

    if (boost::in_degree(u, g) > 0) {
      double probContext = std::accumulate(m_condProbs.begin(), m_condProbs.end(), 1.0, std::multiplies<double>());
      m_C -= probContext * std::log2(probContext);

#ifndef NDEBUG
      std::cout << "Conditional probabilities: ";
      std::copy(m_condProbs.begin(), m_condProbs.end(), std::ostream_iterator<double>(std::cout, " "));
      std::cout << std::endl;
      std::cout << "Probability of the context: " << probContext << std::endl;
      std::cout << "Intermediate result of C: " << m_C << std::endl;
#endif /* NDEBUG */

      double entropy = 0.0;
      BOOST_FOREACH(Edge e, boost::out_edges(u, g)) {
        double condProb = boost::lexical_cast<double>(m_freqMap[boost::target(e, g)])
            / boost::lexical_cast<double>(m_freqMap[boost::source(e, g)]);
        entropy -= condProb * std::log2(condProb);
        m_bic +=  boost::lexical_cast<double>(m_freqMap[boost::target(e, g)]) * std::log(condProb);
      }
      m_h += probContext * entropy;

#ifndef NDEBUG
      std::cout << "Entropy of conditionals: " << entropy << std::endl;
      std::cout << "Intermediate result of h: " << m_h << std::endl;
      std::cout << "Intermediate result of bic: " << m_bic << std::endl;
#endif /* NDEBUG */
    }
  }

  template <typename Edge, typename Graph>
  void examine_edge(Edge e, Graph &g) {
#ifndef NDEBUG
    std::cout << "Examine edge: " << g[boost::source(e, g)].name << "-" << g[boost::target(e, g)].name << std::endl;
#endif /* NDEBUG */

    double condProb = boost::lexical_cast<double>(m_freqMap[boost::target(e, g)])
        / boost::lexical_cast<double>(m_freqMap[boost::source(e, g)]);

#ifndef NDEBUG
    std::cout << "Target Frequency: " << boost::lexical_cast<double>(m_freqMap[boost::target(e, g)]) << std::endl;
    std::cout << "Source Frequency: " << boost::lexical_cast<double>(m_freqMap[boost::source(e, g)]) << std::endl;
    std::cout << "Conditional Probability: " << condProb << std::endl;
#endif /* NDEBUG */

    m_condProbs.push_back(condProb);
  }

  template <typename Edge, typename Graph>
  void tree_edge(Edge e, Graph &g) {
#ifndef NDEBUG
    std::cout << "Tree edge: " << g[boost::source(e, g)].name << "-" << g[boost::target(e, g)].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Edge, typename Graph>
  void back_edge(Edge e, Graph &g) {
#ifndef NDEBUG
    std::cout << "Back edge: " << g[boost::source(e, g)].name << "-" << g[boost::target(e, g)].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Edge, typename Graph>
  void forward_or_cross_edge(Edge e, Graph &g) {
#ifndef NDEBUG
    std::cout << "Forward or Cross edge: " << g[boost::source(e, g)].name << "-" << g[boost::target(e, g)].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Edge, typename Graph>
  void finish_edge(Edge e, Graph &g) {
#ifndef NDEBUG
    std::cout << "Finish edge: " << g[boost::source(e, g)].name << "-" << g[boost::target(e, g)].name << std::endl;
#endif /* NDEBUG */
  }

  template <typename Vertex, typename Graph>
  void finish_vertex(Vertex u, Graph &g) {
#ifndef NDEBUG
    std::cout << "Finish vertex: " << g[u].name << std::endl;
#endif /* NDEBUG */

    m_condProbs.pop_back();
  }

 private:
  VertexFreqMap &m_freqMap;
  double &m_h;
  double &m_C;
  double &m_bic;
  std::vector<double> m_condProbs;
};


}
}

#endif
