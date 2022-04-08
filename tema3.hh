#ifndef tema3
#define tema3

#define EX1 1

void readFromFile(const char* path, const int rank, std::vector<std::pair<int, std::vector<int>>> &perechi);
void printVectorSorted(const std::vector<std::pair<int, std::vector<int>>> &v, const int &rank);
void sendRank(const int &rankSend, const int &source, const int rank);
void sendSize(const int &vectorSize, const int &source, const int rank);
void sendVector(int* vector, const int &vectorSize, const int &source, const int rank);
void receiveRank(int &rankRecv, const int rank);
void receiveSize(int &vectorSize, const int rank);
std::vector<int> receiveVector(const int &vectorSize, const int rank);
int* c_array(const std::vector<int>& v, const int &l = 0, int r = -1);
std::vector<int> generateVector(const int &size);
void sendTopology(int rank, const int destination, const std::vector<std::pair<int, std::vector<int>>> &perechi);
void receiveTopology(const int &source, std::vector<std::pair<int, std::vector<int>>> &perechi);
void sendRankToWorkers(const int &rank, const std::vector<std::pair<int, std::vector<int>>> &perechi);
void sendTopologyToWorkers(const int &rank, const std::vector<std::pair<int, std::vector<int>>> &perechi);
void assembleVector(const int &vectorDimension, const int &chunk, const std::vector<std::pair<int, std::vector<int>>> &perechi);
void solveRank0(char *argv[], const int &numtasks, const int &rank, std::vector<std::pair<int, std::vector<int>>> perechi);
void solveRank1(char *argv[], int numtasks, int rank, std::vector<std::pair<int, std::vector<int>>> perechi);
void solveRank2(char *argv[], int numtasks, int rank, std::vector<std::pair<int, std::vector<int>>> perechi);
void solveRankX(char *argv[], int numtasks, int rank, std::vector<std::pair<int, std::vector<int>>> perechi);
#endif
