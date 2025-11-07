#pragma once

#include "utils.h"

using namespace std;

namespace RankUpAlgorithm {

  class VecN {
  private:
    string _literal;
    LD _TF;
    LD _IDF;
    LD _TF_IDF;
    
  public:
    VecN(string term, int freq, int total_terms)
    : _TF(LD(freq) / LD(total_terms)), _literal(term) {}

    string literal() const { return this->_literal; } 

    LD TF() const { return this->_TF; } 
    LD IDF() const { return this->_IDF; } 
    LD TF_IDF() const { return this->_TF_IDF; } 

    void set_IDF(LD _IDF) { this->_IDF = _IDF; }
    void set_TF_IDF(LD _TF_IDF) { this->_TF_IDF = _TF_IDF; }
  };

  class Document {
  private:
    vector<VecN> vectors;
    map<string, int> terms_freq;

    string _path;
    string _literal;

    LD _magnitude = 0;

    void calc_TF(map<string, unordered_set<string>>& total_terms_freq) {
      ifstream file(this->_path);

      if (!file.is_open()) {
        return;
      }

      string s;
      while (file >> s) {
        if (s.empty()) continue;

        string tl_s = s;
        transform(tl_s.begin(), tl_s.end(), tl_s.begin(), ::tolower);

        terms_freq[tl_s]++;
        total_terms_freq[tl_s].emplace(this->_path);
        
        this->_literal += s + " ";
      }

      if (this->_literal.size()) {
        this->_literal.pop_back();
      }
      
      for (auto& [term, freq] : terms_freq) {
        VecN vec(term, freq, terms_freq.size());
        vectors.emplace_back(vec);
      }
    }
    
  public:
    Document(string path, map<string, unordered_set<string>>& total_terms_freq) {
      this->_path = path;
      this->calc_TF(total_terms_freq);
    }
    
    Document(string literal, map<string, int>& terms_freq, vector<VecN>& vectors)
    : _literal(literal), terms_freq(terms_freq), vectors(vectors) {}

    Document() {}

    void calc_IDF(int total_documents, map<string, unordered_set<string>>& total_terms_freq) {
      LD TF_acc = 0;
      
      for (auto& vec : this->vectors) {
        auto it = total_terms_freq.find(vec.literal());

        if (it == total_terms_freq.end() || it->second.empty()) continue;

        LD docs_with_term = it->second.size();
        LD IDF = log2l(1.0L + LD(total_documents) / docs_with_term);

        vec.set_IDF(IDF);
        vec.set_TF_IDF(vec.TF() * vec.IDF());

        TF_acc += vec.TF_IDF() * vec.TF_IDF();
      }

      this->_magnitude = sqrtl(TF_acc);
    }

    vector<VecN> get_vectors() const { return this->vectors; }
    string literal() const { return this->_literal; }
    LD magnitude() const { return this->_magnitude; }

  };

  class Engine {
  private:
    vector<Document> documents;
    map<string, unordered_set<string>> total_terms_freq;
    bool _ready = false;
    
  public:
    Engine(string dir_path) {
      this->read_documents(dir_path);

      if (documents.size()) {
        this->calc_documents_IDF();
      }
    }

    void calc_documents_IDF() {
      for (auto& doc : this->documents) {
        doc.calc_IDF(this->documents.size(), total_terms_freq);
      }
    }

    void read_documents(string dir_path) {
      for (const auto& file : filesystem::directory_iterator(dir_path)) {
        Document doc(file.path().string(), this->total_terms_freq);

        if (doc.get_vectors().empty()) continue;

        documents.emplace_back(doc);
      }

      if (this->documents.size()) {
        this->_ready = true;
      }
    }

    vector<pair<LD, Document>> search(string query) {
      map<string, int> query_terms_freq;
      string query_literal;
      stringstream ss(query);
      
      while (ss >> query) {
        string tl_str = query;
        transform(tl_str.begin(), tl_str.end(), tl_str.begin(), ::tolower);

        query_terms_freq[tl_str]++;
        query_literal += query + " ";
      }

      query_literal.pop_back();

      vector<VecN> query_vectors;

      LD query_doc_magnitude = 0;
      
      for (auto& [term, freq] : query_terms_freq) {
        VecN vec(term, freq, query_terms_freq.size());
        query_vectors.emplace_back(vec);
      }

      Document query_doc(query_literal, query_terms_freq, query_vectors);
      query_doc.calc_IDF(this->documents.size(), this->total_terms_freq);

      if (query_doc.magnitude() == 0) {
        return {};
      }
      
      vector<pair<LD, Document>> result;
      
      for (const auto& doc : this->documents) {
        
        LD esc_prod = 0;

        for (const auto& query_vec : query_doc.get_vectors()) {
          
          for (const auto& vec : doc.get_vectors()) {
            
            if (query_vec.literal() == vec.literal()) {

              esc_prod += query_vec.TF_IDF() * vec.TF_IDF();

            }
          }
        }

        LD sim = esc_prod / (doc.magnitude() * query_doc.magnitude());

        if (sim != 0) {
          result.emplace_back(sim, doc);
        }
      }

      sort(result.begin(), result.end(), [] (pair<LD, Document>& p1, pair<LD, Document>& p2) {
        return p1.first > p2.first;
      });

      return result;
    }
    
    bool ready() const noexcept { return this->_ready; }

    void DEBUG() {
      for (const auto& doc : this->documents) {
        println("\nDocument: {}", doc.literal());
        for (const auto& vec : doc.get_vectors()) {
          println("{}: TF={:.2f} | IDF={:.2f} | TF-IDF={:.2f}", vec.literal(), vec.TF(), vec.IDF(), vec.TF_IDF());
        }
      }
      println();
    }

  };
  
}