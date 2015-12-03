/*
 *  Bron & Kerbosch Algorithm to find Maximal K-Plexes
 *
 *  Reference: Algorithm 457
 *  Finding All Cliques of an Undirected Graph[H]
 *  Communications of the ACM, September 1973, Volume 16, Number 9
 *
 *  Created by Cynthia Wood on 6/?/2012
 *  Adapted by Tim Becker-Fall 2015
 *
 */

#include "modcppstd.hh"

#define inf 10000000

class Operations{
  int n;
  
public:
  int size;
  int *comp;
  Operations(int maxSize);
  ~Operations();
  void CompAlloc(int);
  void Add(int);
  void Delete();
  void Remove(int);
  void OutputK_Plex();
  vector<int> Duplicate();
};

Operations::Operations(int maxSize){
  n = maxSize;
  // we allocate the data when the object is being created
  comp = new int[n];
  size = 0;
  
}

Operations::~Operations(){
  delete[] comp;
}

void Operations::Add(int AddSel){
  comp[size++] = AddSel;
  
}

void Operations::Delete(void){
  size--;
}

void Operations::Remove(int i){
  
  int j;	
  int index = 0;
  
  for (j = 0; j < size; j++)
    if(comp[j] != i)
      comp[index++] = comp[j];
  size--;
}

void Operations::OutputK_Plex(){
  cerr<<"size: "<<size<<endl;
  for(int i = 0; i < size; i++)
    if(i < size-1)
      cerr<<"   "<<comp[i];
    else{
      cerr<<"   "<<comp[i]<<endl;
    }
}

vector<int> Operations::Duplicate(){

  vector<int> k_plex;
  for(int i = 0; i < size; i++){
    k_plex.push_back(comp[i]);
  }

  return k_plex;


}

class Graph{
  int n;
  int **C;
  int* degrees;
  Operations *compsub;
 
public:
  Graph(char *InputFile, int kk);
  // ~Graph();
  void MaxK_Plexes(int kk);
  int size(void);
  void OutputMatrix(void);
  void PrintDot(int prefix, int *shownNodes, int shownNodesCount);
  int Update(int start, int end, int kk, int *old, int nc, int *New, int s);
private:
  void extension(int *old, int ne, int ce, int kk);
  bool inSet(int node, int *set, int setSize);
  
};

Graph::Graph(char *InputFile, int kk){

  ifstream f_in;
  f_in.open(InputFile);
  int m;

  f_in>>n;
  f_in>>m;

  int *x = new int[m];
  int *y = new int[m];
  int z;

  for(int i = 0; i < m; i++){
    f_in>>x[i]>>y[i]>>z;
  }


  //Create matrix
  C = new int*[n];
  for( int i = 0; i < n;i++)
    C[i] = new int[n];
  
  for(int i = 0; i < n-1; i++){
    for(int j = i; j < n; j++) C[i][j] = 0;
  }

  for(int i = 0; i < m; i++){
    C[x[i]][y[i]] = 1;
    C[y[i]][x[i]] = 1;
  }

  for(int i = 0; i < n; i++){
    C[i][i] = 1;
  }
 
  compsub = new Operations(n);
  f_in.close();

}

/*Graph::~Graph(){
  
  for(int i = 0; i < n; i++)
    delete[] C[i];
  delete[] C;
  delete compsub;
  delete Neighbors;
  }*/

void Graph::MaxK_Plexes(int kk){

  int *all = new int[n];
  for(int i = 0; i < n; i++){
    all[i] = i;
  }

  extension(all, 0, n,kk);
  delete[] all;
}

int Graph::size(){
  return n;
}

bool Graph::inSet(int node, int set[], int setSize) {
  for (int i = 0; i < setSize; i++) {
    if (set[i] == node) {
      return 1;
    }
  }
  return 0;
}

void Graph::PrintDot(int prefix, int shownNodes[], int shownNodesCount) {
  cout << "subgraph A_";
  cout << prefix;
  cout << " {\n";
  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++){
      if(C[i][j] && i != j && inSet(i, shownNodes, shownNodesCount)) {
	cout << prefix;
	cout << "_";
	cout << i;
	cout << " -- ";
	cout << prefix;
	cout << "_";
	cout << j;
	cout << "\n";
      }
    }
  }
  cout << "}\n";
  
}

int Graph::Update(int start, int end, int kk, int old[], int nc, int New[], int s){


  int p;
  vector<int> k_plex = compsub->Duplicate();
  int *degrees = new int[k_plex.size()];

  for(int i = 0; i < k_plex.size(); i++){
    degrees[i] = 0;
    for(int j = 0; j < k_plex.size(); j++){
      if(C[k_plex.at(i)][k_plex.at(j)] == 1 && i != j){
	degrees[i]++;
      }
    }
  }

  //FIND NEXT NODE FOR CANDIDATES/NOT
  for(int i = start; i < end; i++){
    p = old[i];
    for(int j = 0; j < k_plex.size(); j++){
      if(C[k_plex.at(j)][p] == 1){
	degrees[j]++;
      }
    }
    int min = inf;
    for(int q = 0; q < k_plex.size(); q++){
      if(degrees[q] < min){
	min = degrees[q];
      }
    }
    


    //CHECK IF NEWEST NODE KEEPS K-PLEX OR NOT
    if(min >= int(k_plex.size() + 1 - kk) && p != s){
      New[nc] = old[i];
      nc++;
    }
    
    for(int j = 0; j < k_plex.size(); j++){
      if(C[k_plex.at(j)][p] == 1){
	degrees[j]--;
      }
    }
  }

  delete[] degrees;
  return nc;
}


//Extension operator to find all the configurations of the compsub
void Graph::extension(int old[], int ne, int candidates, int kk){
  
  //Bactrackcycle
  int *New = new int[candidates];
  int s, newne, newCandidates;
  //CHANGE TO UNION OF NEIGHBORS OF CURRENT K-PLEX
  for(int i = ne; i < candidates; i++){
    newne = 0;
    s = old[i];
    compsub->Add(s);

    //Fill new set not
    newne = Update(0,ne,kk, old,newne,New,s);
    
    newCandidates = newne;
    //Fill new set candidates
    newCandidates = Update(ne,candidates,kk,old,newCandidates,New,s);

   
    if(newCandidates == 0){ 
      compsub->OutputK_Plex();
    }
    else if(newne < newCandidates){
       extension(New, newne, newCandidates,kk);
    }
    
     //Remove from compsub
    compsub->Delete();
    
    //Add to not
    ne++;
     //End Backtrackcycle
  }
  
  delete[] New;		
}

int main(int argc, char *argv[]){
  
  ifstream f_in;
  f_in.open(argv[1]);
  if(!f_in.is_open()){
    cout << "CANNOT OPEN FILE" << endl;
  }
  else{
    Graph g(argv[1],atoi(argv[2]));
    clock_t begin = clock();
    g.MaxK_Plexes(atoi(argv[2]));
    clock_t end = clock();

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout <<"Time to solve: "<<elapsed_secs<<"s";
    double elapsed_mins = elapsed_secs/60;
    cout << "(or " << elapsed_mins<<"m)"<<endl;
  }


  return 0;
  }
