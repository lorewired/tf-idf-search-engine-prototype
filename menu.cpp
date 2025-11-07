#include "engine.hpp"

int main() {
  using namespace std;
  using namespace RankUpAlgorithm;

  print("Digite o diretório para alimentação: ");
  string dir;
  cin >> dir;

  RankUpAlgorithm::Engine engine(dir);

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

    if (query[0] == '0') {
      break;
    }
    
    vector<pair<LD, Document>> result = engine.search(query);

    if (result.empty()) {
      println("\nResultado: Sua pesquisa não encontrou nenhum documento correspondente.");
      print("\nPressione qualquer tecla para continuar...");
      getchar();
      continue;
    }

    for (int i = 0; i < result.size(); i++) {
      println("\nResultado {} [SCORE {:.2f}] : {}", i + 1, result[i].first, result[i].second.literal());
    }
    
    print("\nPressione qualquer tecla para continuar...");

    getchar();
  }
  
  return 0;
}