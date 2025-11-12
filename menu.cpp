#include "engine.hpp"

int main() {
  using namespace std;
  using namespace rank_up_algorithm;

  rank_up_algorithm::engine engine("feeding");

  if (!engine.ready()) {
    system("clear");
    println("Falha ao ler os documentos de alimentação.");
    return 0;
  }

  while (true) {
    system("clear");
    
    print("[0 para sair] Pesquisar: ");

    string query;
    while (getline(cin, query) && query.empty());

    if (query[0] == '0') break;
    
    pair<vec, vector<pair<LD, vec>>> results = engine.search(query);

    if (results.second.empty()) {
      println("\nResultado: Sua pesquisa não encontrou nenhum documento correspondente.");
      print("\nPressione qualquer tecla para continuar...");
      getchar();
      continue;
    }

    for (int i = 0; i < results.second..size(); i++) {
      println("\nResultado {} [SCORE {:.2f}] : {}", i + 1, results.second[i].first, results.second[i].second.content());
    }

    print("\nPressione qualquer tecla para continuar...");

    getchar();
  }
  
  return 0;
}
