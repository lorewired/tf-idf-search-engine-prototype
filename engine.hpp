#pragma once

#include "utils.h"
#include <filesystem>

using namespace std;

namespace rank_up_algorithm {
  constexpr int SUFFIX_ERASE_LIMIT = 2;

  class vec_term {
    string _literal;
    LD _TF = 0.0l;
    LD _IDF = 0.0l;
    LD _TF_IDF = 0.0l;

  public:
    vec_term(string __literal, int freq, int total_terms)
    : _TF(LD(freq) / LD(total_terms)), _literal(__literal) {}

    vec_term() {}

    string literal() const { return this->_literal; }
    LD TF() const { return this->_TF; } 
    LD IDF() const { return this->_IDF; } 
    LD TF_IDF() const { return this->_TF_IDF; } 

    void set_IDF(LD _IDF) { this->_IDF = _IDF; }
    void set_TF_IDF(LD _TF_IDF) { this->_TF_IDF = _TF_IDF; }

  };
  
  class pref_tree {
    class pref_tree_node {
      vec_term _term;
      
    public:
      pref_tree_node() {}

      vector<vec_term> _suffixes;
      map<char, pref_tree_node *> _adjacents;

      vec_term term() const { return this->_term; }

      void set_term(const vec_term& __term) { this->_term = __term; }
    };

    pref_tree_node * _root = nullptr;

  public:
    pref_tree() : _root(new pref_tree_node) {}

    void add(const vec_term& new_term, string cur_str = "", pref_tree_node * cur = nullptr) {
      if (cur == nullptr) cur = this->_root;

      if (cur_str == new_term.literal()) {
        cur->_suffixes.push_back(new_term);
        cur->set_term(new_term);
        return;
      }
      
      char target = new_term.literal()[cur_str.size()];

      auto it = cur->_adjacents.find(target);

      pref_tree_node * next = it == cur->_adjacents.end() ? new pref_tree_node : it->second;
      next->set_term(new_term);

      cur->_adjacents[target] = next;

      add(new_term, cur_str + target, next);
      
      if (cur != this->_root) {
        cur->_suffixes.push_back(new_term);
      }
    }

    vector<vec_term> matches(const string& target_str) const {
      pref_tree_node* cur = this->_root;
      string cur_str;

      stack<pref_tree_node *> stk;

      for (char c : target_str) {
        auto it = cur->_adjacents.find(c);
        if (it == cur->_adjacents.end()) break;

        stk.push(it->second);
        
        cur_str += c;
        cur = it->second;
      }

      unordered_map<string, vec_term> mp_matches;

      auto add_matches = [&] (pref_tree_node * node, LD additional_weight) {
        for (const vec_term& suff : node->_suffixes) {
          if (mp_matches.find(suff.literal()) != mp_matches.end()) {
            continue;
          }

          vec_term tmp_suff = suff;

          int dist = abs(int(target_str.size() - tmp_suff.literal().size()));

          if (dist > SUFFIX_ERASE_LIMIT) continue;

          LD penalty = 5.0L / (5.0L + dist + additional_weight);

          tmp_suff.set_TF_IDF(tmp_suff.TF_IDF() * penalty);

          auto& cur_match = mp_matches[tmp_suff.literal()];
          if (cur_match.TF_IDF() == 0 || tmp_suff.TF_IDF() > cur_match.TF_IDF()) {
            cur_match = tmp_suff;
          }
        }
      };

      if (cur_str == target_str) {
        add_matches(cur, 0);
        stk.pop();
      }

      LD additional_weight = 0;

      if (stk.size()) {
        int initial_dist = target_str.size() - cur_str.size();
        if (initial_dist > SUFFIX_ERASE_LIMIT) return {};
        
        additional_weight = LD(initial_dist);
      }

      for (int i = 0; i < SUFFIX_ERASE_LIMIT && stk.size(); i++) {
        pref_tree_node * prev = stk.top();
        stk.pop();

        add_matches(prev, additional_weight);

        additional_weight += 1.0L;
      }

      vector<vec_term> all_matches;

      for (const auto& [_, match] : mp_matches) {
        all_matches.push_back(match);
      }

      return all_matches;
    }
  };

  class vec {
    pref_tree _prefixes;
    vector<vec_term> _terms;
    unordered_map<string, int> _terms_freq;
    
    LD _magnitude = 0;
    string _content;

  public:
    vec() {}
    
    vec(
      const string& file_path,
      unordered_map<string, unordered_set<string>>& corpus_terms_freq
    ) {
      stringstream ss;

      bool is_file = filesystem::exists(file_path);
      
      if (is_file) {
        ifstream file(file_path);

        if (!file.is_open()) return;

        ss << file.rdbuf();
      } else {
        string content = file_path;
        ss = stringstream(content);
      }
      
      string term_str;

      while (ss >> term_str) {
        if (term_str.empty()) continue;

        string lower_case_term = term_str;
        transform(lower_case_term.begin(), lower_case_term.end(), lower_case_term.begin(), ::tolower);

        if (is_file) {
          corpus_terms_freq[lower_case_term].insert(file_path);
        }

        this->_terms_freq[lower_case_term]++;
        this->_content += term_str + ' ';
      }

      if (this->_content.size()) this->_content.pop_back();

      for (const auto& [term_str, freq] : this->_terms_freq) {
        vec_term term(term_str, freq, this->_terms_freq.size());
        this->_terms.push_back(term);
      }
    }

    void calc_IDF(const int& vectors_total, const unordered_map<string, unordered_set<string>>& corpus_terms_freq) {
      LD TF_IDF_acc = 0;

      for (vec_term& term : this->_terms) {
        auto it = corpus_terms_freq.find(term.literal());

        LD vectors_with_term = 0;

        if (it != corpus_terms_freq.end()) {
          vectors_with_term = it->second.size();
        }

        LD IDF = log2l(1.0L + LD(vectors_total) / (vectors_with_term + 1.0L));

        term.set_IDF(IDF);
        term.set_TF_IDF(term.TF() * term.IDF());

        TF_IDF_acc += term.TF_IDF() * term.TF_IDF();

        this->prefixes().add(term);
      }

      this->_magnitude = sqrtl(TF_IDF_acc);
    }

    pref_tree& prefixes() { return this->_prefixes; }
    const pref_tree& prefixes() const { return this->_prefixes; }
    vector<vec_term> terms() const { return this->_terms; }
    unordered_map<string, int> terms_freq() const { return this->_terms_freq; }
    LD magnitude() const { return this->_magnitude; }
    string content() const { return this->_content; }

  };

  class corpus {
    unordered_map<string, unordered_set<string>> _corpus_terms_freq;
    vector<vec> _vectors;

  public:
    bool read_vectors(const string& __dir_path) {
      for (const auto& file : filesystem::directory_iterator(__dir_path)) {
        vec doc(file.path().string(), this->_corpus_terms_freq);

        if (doc.content().size()) {
          this->_vectors.push_back(doc);
        }
      }
      return this->_vectors.size();
    }

    void calc_terms_IDF() {
      for (auto& vec : this->_vectors) {
        vec.calc_IDF(this->_vectors.size(), this->_corpus_terms_freq);
      }
    }

    vector<vec> vectors() const { return this->_vectors; }
    unordered_map<string, unordered_set<string>> corpus_terms_freq() const { return this->_corpus_terms_freq; }
  };

  class engine {
    corpus corp;
    bool _ready = false;
  public:
    engine(const string& __dir_path) {
      this->_ready = this->corp.read_vectors(__dir_path);

      if (this->_ready) {
        this->corp.calc_terms_IDF();
      }
    }

    pair<vec, vector<pair<LD, vec>>> search(const string& query_str) {
      if (query_str.empty()) return make_pair(vec(), vector<pair<LD, vec>>());

      unordered_map<string, unordered_set<string>> ignore;
      
      vec query_vec(query_str, ignore);
      query_vec.calc_IDF(this->corp.vectors().size(), this->corp.corpus_terms_freq());

      vector<pair<LD, vec>> result;

      for (const auto& vec : this->corp.vectors()) {
        if (vec.magnitude() == 0) continue;

        LD esc_prod = 0;

        for (const vec_term& query_term : query_vec.terms()) {
          vector<vec_term> vec_prefixes = vec.prefixes().matches(query_term.literal());

          if (vec_prefixes.empty()) continue;

          vec_term best_match = *max_element(
            vec_prefixes.begin(),
            vec_prefixes.end(),
            [] (const vec_term& term1, const vec_term& term2) { return term1.TF_IDF() < term2.TF_IDF(); }
          );

          esc_prod += query_term.TF_IDF() * best_match.TF_IDF();
        }

        LD sim = esc_prod / (vec.magnitude() * query_vec.magnitude());

        if (sim != 0) {
          result.emplace_back(sim, vec);
        }
      }
      
      sort(result.begin(), result.end(), [] (const pair<LD, vec>& p1, const pair<LD, vec>& p2) {
        return p1.first > p2.first;
      });

      return {query_vec, result};
    }

    void DEBUG() {
      for (const vec& vec : this->corp.vectors()) {
        println("\nDocument: {}", vec.content());
        for (const auto& term : vec.terms()) {
          println("{}: TF={:.2f} | IDF={:.2f} | TF-IDF={:.2f}", term.literal(), term.TF(), term.IDF(), term.TF_IDF());
        }
      }
    }

    bool ready() const { return this->_ready; }
  };
}
