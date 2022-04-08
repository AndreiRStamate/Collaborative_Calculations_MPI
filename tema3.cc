#include "mpi.h"
#include <stdlib.h>
#include <cstdio>
#include <vector>
#include "tema3.hh"
#include <algorithm>
#include <sstream>
#include <iostream>


int main(int argc, char *argv[]) {
    std::vector<std::pair<int, std::vector<int>>> perechi;

    int numtasks, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   

    if (rank == 0) {
        solveRank0(argv, numtasks, rank, perechi);

    } else if (rank == 1) {
        solveRank1(argv, numtasks, rank, perechi);

    } else if (rank == 2){
        solveRank2(argv, numtasks, rank, perechi);

    } else {
        solveRankX(argv, numtasks, rank, perechi);
    }

    MPI_Finalize();
}


void readFromFile(const char* path, const int rank, std::vector<std::pair<int, std::vector<int>>> &perechi) {
    auto fp = freopen(path, "r", stdin);

    int workerCount;
    auto _ = scanf("%d", &workerCount);
    std::vector<int> workers(workerCount);
    for (int i = 0; i < workerCount; i++) {
        _ = scanf("%d", &workers[i]);
    }
    perechi.emplace_back(std::make_pair(rank, workers));

    fclose(fp);
}

void printVectorSorted(const std::vector<std::pair<int, std::vector<int>>> &v, const int &rank) {
    std::vector<std::pair<int, std::vector<int>>> perechiSortat(3);
    // aranjare perechi dupa rang
    for (int i = 0; i < v.size(); i++) {
        perechiSortat[v[i].first] = std::make_pair(v[i].first, v[i].second);
    }
    // afisare conform cerintei
    std::stringstream s;
    s << rank << " -> ";
    for (auto pereche : perechiSortat) {
        s << pereche.first << ":";
        int cnt = 0;
        for (auto e : pereche.second) {
            s << e;
            if (cnt++ < pereche.second.size()-1) s << ","; 
        }
        s << " ";
    }
    s << std::endl;
    std::cout << s.rdbuf();
}

void sendRank(const int &rankSend, const int &source, const int rank) {
    // se trimite rank-ul catre procesul rank
    if (EX1) std::cout << "M(" << source << "," << rank << ")" << std::endl;
    MPI_Send(&rankSend, 1, MPI_INT, rank, 69, MPI_COMM_WORLD);
}

void sendSize(const int &vectorSize, const int &source, const int rank) {
    // se trimite marimea catre procesul rank
    sendRank(vectorSize, source, rank);
}

void sendVector(int* vector, const int &vectorSize, const int &source, const int rank) {
    // se trimite un vector catre procesul rank
    if (EX1) std::cout << "M(" << source << "," << rank << ")" << std::endl; 
    MPI_Send(vector, vectorSize, MPI_INT, rank, 69, MPI_COMM_WORLD);
}

void receiveRank(int &rankRecv, const int rank) {
    // se primeste rankul de la procesul rank
    MPI_Recv(&rankRecv, 1, MPI_INT, rank, 69, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void receiveSize(int &vectorSize, const int rank) {
    // se primeste marimea de la procesul rank
    receiveRank(vectorSize, rank);
}

std::vector<int> receiveVector(const int &vectorSize, const int rank) {
    int* vectorRecv = (int*) malloc(vectorSize * sizeof(int));

    // se primeste vectorul de la procesul rank
    MPI_Recv(vectorRecv, vectorSize, MPI_INT, rank, 69, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
    
    std::vector<int> v(vectorRecv, vectorRecv + vectorSize);
    free(vectorRecv);

    return v;
}

int* c_array(const std::vector<int>& v, const int &l, int r) {
    if (r == -1) r = v.size();
    int* arr = (int*) malloc((r-l) * sizeof(int));
    int k = 0;
    for (int i = l; i < r; i++) {
        arr[k++] = v[i];
    }
    return arr;
}

std::vector<int> generateVector(const int &size) {
    std::vector<int> v(size);
    for (int i = 0; i < size; i++) {
        v[i] = i;
    }
    return v;
}

void sendTopology(int rank, const int destination, const std::vector<std::pair<int, std::vector<int>>> &perechi) {
    // trimitere dimensiuni cluster0
    int vectorSize = perechi[0].second.size();
    sendSize(vectorSize, rank, destination);

    // trimitere topologie catre coordonatori
    int* vectorToPointer = c_array(perechi[0].second);
    sendVector(vectorToPointer, vectorSize, rank, destination);
    free(vectorToPointer);
}

void sendTopologyBonus(int rank, const int destination, const std::vector<std::pair<int, std::vector<int>>> &perechi) {
    // trimitere dimensiuni cluster0
    int index = destination == 0 ? 2 : 1;
    int vectorSize = perechi[index].second.size();
    sendSize(vectorSize, rank, destination);

    // trimitere topologie catre coordonatori
    int* vectorToPointer = c_array(perechi[index].second);
    sendVector(vectorToPointer, vectorSize, rank, destination);
    free(vectorToPointer);
}

void receiveTopology(const int &source, std::vector<std::pair<int, std::vector<int>>> &perechi) {
    //primire dimensiuni + date
    int vectorSize;
    receiveSize(vectorSize, source);
    perechi.emplace_back(std::make_pair(source, receiveVector(vectorSize, source)));
}

void receiveTopologyBonus(const int &source, const int index, std::vector<std::pair<int, std::vector<int>>> &perechi) {
    //primire dimensiuni + date
    int vectorSize;
    receiveSize(vectorSize, source);
    perechi.emplace_back(std::make_pair(index, receiveVector(vectorSize, source)));
}

void sendRankToWorkers(const int &rank, const std::vector<std::pair<int, std::vector<int>>> &perechi) {
    // trimitere rank catre workeri
    for (int i = 0; i < perechi[0].second.size(); i++) {
        sendRank(rank, rank, perechi[0].second[i]);
    }
}

void sendTopologyToWorkers(const int &rank, const std::vector<std::pair<int, std::vector<int>>> &perechi) {
    // trimitere marimi vectori topologie catre wokeri
    for (int i = 0; i < perechi[0].second.size(); i++) {
        sendSize(perechi.size(), rank, perechi[0].second[i]);
    }

    // trimitere topologie catre workeri
    for (int i = 0; i < perechi.size(); i++) {
        for (int j = 0; j < perechi[0].second.size(); j++) {
            sendRank(perechi[i].first, rank, perechi[0].second[j]);
            int auxSize = perechi[i].second.size();
            sendSize(auxSize, rank, perechi[0].second[j]);
            sendVector(c_array(perechi[i].second), auxSize, rank, perechi[0].second[j]);
        }
    }
}

void assembleVector(const int &vectorDimension, const int &chunk, const std::vector<std::pair<int, std::vector<int>>> &perechi) {
    std::vector<int> finalVector(vectorDimension);
    int cnt = 0;

    // bucata cluster0
    for (int i = 0; i < perechi[0].second.size(); i++) {
        std::vector<int> auxVector0 = receiveVector(chunk, perechi[0].second[i]);
        for (auto j : auxVector0) {
            finalVector[cnt++] = j;
        }
    }

    // bucata cluster1
    std::vector<int> auxVector1 = receiveVector(vectorDimension, 2);
    for (int i = chunk * perechi[0].second.size(); 
        i < chunk * (perechi[0].second.size() + perechi[1].second.size()); i++) {
        finalVector[i] = auxVector1[i];
    }

    // bucata cluster2
    std::vector<int> auxVector2 = receiveVector(vectorDimension, 2);
    for (int i = chunk * (perechi[0].second.size() + perechi[1].second.size()); 
        i < vectorDimension; i++) {
        finalVector[i] = auxVector2[i];
    }

    std::stringstream s;
    s << "Rezultat: ";
    for (int i = 0; i < finalVector.size(); i++) s << finalVector[i] << " ";
    s << std::endl;

    std::cout << s.rdbuf();
}

void solveRank0(char *argv[], const int &numtasks, const int &rank, std::vector<std::pair<int, std::vector<int>>> perechi) {
    int vectorDimension = atoi(argv[1]);

    // citire date fisier
    readFromFile("cluster0.txt", 0, perechi);

    //sendTopology(rank, 1, perechi);
    sendTopology(rank, 2, perechi);

    //receiveTopology(1, perechi);
    receiveTopology(2, perechi);

    receiveTopologyBonus(2, 1, perechi);

    std::vector<std::pair<int, std::vector<int>>> perechiBonus(3);
    perechiBonus[0].first = perechi[0].first;
    for (int i = 0; i < perechi[0].second.size(); i++) {
        perechiBonus[0].second.push_back(perechi[0].second[i]);
    }

    perechiBonus[1].first = perechi[2].first;
    for (int i = 0; i < perechi[2].second.size(); i++) {
        perechiBonus[1].second.push_back(perechi[2].second[i]);
    }

    perechiBonus[2].first = perechi[1].first;
    for (int i = 0; i < perechi[1].second.size(); i++) {
        perechiBonus[2].second.push_back(perechi[1].second[i]);
    }

    printVectorSorted(perechiBonus, rank);

    sendRankToWorkers(rank, perechiBonus);

    sendTopologyToWorkers(rank, perechiBonus);

    // generare vector de calculat
    std::vector<int> initialVector = generateVector(atoi(argv[1]));
    int initialVectorSize = initialVector.size();

    // trimitere vector de calculat catre procesele coordonator
    sendSize(initialVectorSize, rank, 2);

    int* initialVectorToPointer = c_array(initialVector);
    sendVector(initialVectorToPointer, initialVectorSize, rank, 2);
    free(initialVectorToPointer);

    // trimitere bucati vector catre workeri
    int numberProcesses = perechiBonus[0].second.size() + perechiBonus[1].second.size() + perechiBonus[2].second.size();
    int chunk = vectorDimension / numberProcesses;
    for (int i = 0; i < perechiBonus[0].second.size(); i++) {
        sendSize(chunk, 0, perechiBonus[0].second[i]);
        int* aux = c_array(initialVector, i*chunk, (i+1)*chunk);
        sendVector(aux, chunk, rank, perechiBonus[0].second[i]);
        free(aux);
    }

    // construireVectorDublu
    assembleVector(vectorDimension, chunk, perechiBonus);
}

void solveRank1(char *argv[], int numtasks, int rank, std::vector<std::pair<int, std::vector<int>>> perechi) {
    int vectorDimension = atoi(argv[1]);

    // citire date fisier
    readFromFile("cluster1.txt", 1, perechi);

    //sendTopology(rank, 0, perechi);
    sendTopology(rank, 2, perechi);

    //receiveTopology(0, perechi);
    receiveTopology(2, perechi);

    receiveTopologyBonus(2, 0, perechi);

    std::vector<std::pair<int, std::vector<int>>> perechiBonus(3);
    perechiBonus[0].first = perechi[0].first;
    for (int i = 0; i < perechi[0].second.size(); i++) {
        perechiBonus[0].second.push_back(perechi[0].second[i]);
    }

    perechiBonus[1].first = perechi[2].first;
    for (int i = 0; i < perechi[2].second.size(); i++) {
        perechiBonus[1].second.push_back(perechi[2].second[i]);
    }

    perechiBonus[2].first = perechi[1].first;
    for (int i = 0; i < perechi[1].second.size(); i++) {
        perechiBonus[2].second.push_back(perechi[1].second[i]);
    }

    printVectorSorted(perechiBonus, rank);

    sendRankToWorkers(rank, perechiBonus);

    sendTopologyToWorkers(rank, perechiBonus);

    //primire vector de calculat
    int initialVectorSize;
    receiveSize(initialVectorSize, 2);

    std::vector<int> initialVector = receiveVector(initialVectorSize, 2);

    // trimitere bucati vector catre workeri
    int numberProcesses = perechiBonus[0].second.size() + perechiBonus[1].second.size() + perechiBonus[2].second.size();
    int chunk = vectorDimension / numberProcesses;

    // trimitere la toti workerii
    for (int i = 0; i < perechiBonus[0].second.size(); i++) {
        sendSize(chunk, rank, perechiBonus[0].second[i]);
        int* aux = c_array(initialVector, perechiBonus[1].second.size()*chunk + i*chunk, perechiBonus[1].second.size()*chunk +(i+1)*chunk);
        sendVector(aux, chunk, rank, perechiBonus[0].second[i]);
        free(aux);
    }

    std::vector<int> finalVector(vectorDimension);
    int cnt = perechiBonus[1].second.size()*chunk;
    for (int i = 0; i < perechiBonus[0].second.size(); i++) {
        auto auxVector = receiveVector(chunk, perechiBonus[0].second[i]);
        for (auto j : auxVector) {
            finalVector[cnt++] = j;
        }
    }

    // trimitere vector calculat catre procesul 2
    sendVector(c_array(finalVector), finalVector.size(), rank, 2);
}

void solveRank2(char *argv[], int numtasks, int rank, std::vector<std::pair<int, std::vector<int>>> perechi) {
    int vectorDimension = atoi(argv[1]);

    readFromFile("cluster2.txt", 2, perechi);

    sendTopology(rank, 0, perechi);
    sendTopology(rank, 1, perechi);

    receiveTopology(0, perechi);
    receiveTopology(1, perechi);

    printVectorSorted(perechi, rank);

    sendTopologyBonus(rank, 1, perechi);
    sendTopologyBonus(rank, 0, perechi);

    sendRankToWorkers(rank, perechi);

    sendTopologyToWorkers(rank, perechi);

    //primire vector de calculat
    int initialVectorSize;
    receiveSize(initialVectorSize, 0);

    std::vector<int> initialVector = receiveVector(initialVectorSize, 0);

    // trimitere vector de calculat catre procesele coordonator
    sendSize(initialVectorSize, rank, 1);

    int* initialVectorToPointer = c_array(initialVector);
    sendVector(initialVectorToPointer, initialVectorSize, rank, 1);
    free(initialVectorToPointer);

    //trimitere bucati vector catre workeri
    int numberProcesses = perechi[0].second.size() + perechi[1].second.size() + perechi[2].second.size();
    int chunk = vectorDimension / numberProcesses;
    int chunksSent = chunk * (perechi[1].second.size() + perechi[2].second.size());
    // trimitere la nr_workeri-1
    for (int i = 0; i < perechi[0].second.size()-1; i++) {
        chunksSent += chunk;
        sendSize(chunk, rank, perechi[0].second[i]);
        int* aux = c_array(initialVector, (perechi[1].second.size() + perechi[2].second.size())*chunk + i*chunk, vectorDimension);
        sendVector(aux, chunk, rank, perechi[0].second[i]);
        free(aux);
    }

    //trimitere restul vectorului la ultimul worker
    int remainingChunk = vectorDimension - chunksSent;

    int i = perechi[0].second.size()-1;
    sendSize(remainingChunk, rank, perechi[0].second[i]);
    int* aux = c_array(initialVector, chunksSent, vectorDimension);
    sendVector(aux, remainingChunk, rank, perechi[0].second[i]);
    free(aux);

    std::vector<int> finalVector(vectorDimension);
    // primire nr_workers-1 parti
    int cnt = (perechi[1].second.size() + perechi[2].second.size()) * chunk;
    for (int i = 0; i < perechi[0].second.size()-1; i++) {
        std::vector<int> auxVector = receiveVector(chunk, perechi[0].second[i]);
        for (auto j : auxVector) {
            finalVector[cnt++] = j;
        }
    }
    // primire ultima parte
    std::vector<int> auxVector = receiveVector(remainingChunk, perechi[0].second[i]);
    for (auto j : auxVector) {
        finalVector[cnt++] = j;
    }

    std::vector<int> finalVector2(vectorDimension);
    // bucata cluster1
    std::vector<int> auxVector1 = receiveVector(vectorDimension, 1);
    for (int i = chunk * perechi[1].second.size(); 
        i < chunk * (perechi[1].second.size() + perechi[2].second.size()); i++) {
        finalVector2[i] = auxVector1[i];
    }

    // trimitere vector calculat la procesul 0
    sendVector(c_array(finalVector2), finalVector2.size(), rank, 0);
    sendVector(c_array(finalVector), finalVector.size(), rank, 0);
}

void solveRankX(char *argv[], int numtasks, int rank, std::vector<std::pair<int, std::vector<int>>> perechi) {
    //primire rank coordonator
        int coordRank;
        receiveRank(coordRank, MPI_ANY_SOURCE);

        // primire topologie
        int topologieSize;
        receiveSize(topologieSize, coordRank);
        for (int i = 0; i < topologieSize; i++) {
            int rank;
            int vectorSize;
            receiveRank(rank, coordRank);
            receiveSize(vectorSize, coordRank);
            perechi.emplace_back(std::make_pair(rank, receiveVector(vectorSize, coordRank)));
        }

        // afisare topologie
        printVectorSorted(perechi, rank);

        // primire vector de dublat
        int vectorSize;
        receiveSize(vectorSize, coordRank);

        std::vector<int> aux = receiveVector(vectorSize, coordRank);
        
        // dublare vector
        std::transform(aux.begin(), aux.end(), aux.begin(), [](const int &i) -> int {return i*2;});

        // trimitere vector dublat la procesul coordonator
        sendVector(c_array(aux), aux.size(), rank, coordRank);
}