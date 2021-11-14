/*
    Copyright (C) Fabio Ciani, 2021
    Written by Fabio Ciani <fabio1.ciani@mail.polimi.it>, July 2021
    
    This file is part of the final projects for the B.Sc. in Computer Science & Engineering at Polytechnic University of Milan.
    
    The following piece of software refers to the Algorithms & Data Structures course.
    Its specifications have been presented during the academic year 2020-2021.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

/*
    GLOBAL VARIABLES & STRUCTS
*/

unsigned int d, k;

// GRAPH ELEMENTS
struct Vertex {
    unsigned long long int distance;    // minimum path length to reach node #0
    struct Vertex* previous;    // previous vertex on MST path
    short int isConnectedToSource;
};

typedef struct Vertex vertex;

typedef struct {
    unsigned long long int label;
    vertex* V;      // vertices (i.e., nodes) of the graph
    unsigned int** E;       // (weighted) edges of the graph stored as an adjacency matrix
    unsigned long long int pathsLength;     // total path length (i.e., sum of all minimum path lengths) from node #0 to all the other nodes
} graph;

// MAX-HEAP ELEMENTS
struct Cell {
    unsigned long long int tag, length;     // tag := label, length := pathsLength
};

typedef struct Cell cell;

typedef struct {
    unsigned int heapSize;
    cell* array;    // the stored data will be the actual leaderboard, thus the indexes of the best graphs
} heap;     // max-heap

// PRIORITY-QUEUE ELEMENTS
typedef struct {
    unsigned int heapSize;
    vertex* graphVertices;      // reference to G.V
    unsigned int* vertexIndex;      // the stored data will be the label of a specific (numbered) graph node
    unsigned int* indexPositions;   // contains the index of vertexIndex associated to each (numbered) graph node, hence indexPositions[i] is the index of vertexIndex referred to node i
} queue;    // min-heap

/*
    FUNCTION PROTOTYPES
*/

unsigned int stringToNumber(char*);

void initializeVertex(vertex*);
void initializeGraph(graph*);

void addGraph(graph*);
void freeGraph(graph*);

unsigned int parent(unsigned int);
unsigned int left(unsigned int);
unsigned int right(unsigned int);

void initializeMaxHeap(heap*);

void swapCells(cell*, cell*);

void maxHeapify(heap*, unsigned int);
void heapExtractMax(heap*);
void maxHeapRaiseUp(heap*, unsigned int);
void maxHeapInsert(heap*, cell);
void freeHeap(heap*);

void updateLeaderboard(heap*, graph*);

void initializeQueue(queue*);
void overwriteQueue(queue*, graph*);

void swap(unsigned int*, unsigned int*);

void minHeapify(queue*, unsigned int);
unsigned int heapExtractMin(queue*);
void minHeapInsert(queue*, unsigned int);
void minHeapUpdateKey(queue*, unsigned int);
void freeQueue(queue*);

void optimizedDijkstra(graph* g, queue* q);

void showBestGraphs(heap);

void parseInput();

/*
    FUNCTION DECLARATIONS
*/

unsigned int stringToNumber(char* str) {
    unsigned int n = 0;

    while(*str != '\0') {
        n = ((n << 1) + (n << 3)) + (*str - '0');   // n * 10 = (n * 2) + (n * 8) = (n << 1) + (n << 3)
        str++;      // move pointer to the next character
    }

    return n;
}

void initializeVertex(vertex* x) {
    x->distance = ULLONG_MAX;
    x->previous = NULL;
    x->isConnectedToSource = 0;
}

void initializeGraph(graph* x) {
    unsigned int i;

    vertex* nodes = (vertex *) calloc(d, sizeof(vertex));

    if(nodes == NULL)
        exit(EXIT_FAILURE);

    for(i = 0; i < d; i++)
        initializeVertex(&nodes[i]);

    unsigned int** edges = (unsigned int **) calloc(d, sizeof(unsigned int *));

    if(edges == NULL)
        exit(EXIT_FAILURE);

    for(i = 0; i < d; i++) {
        edges[i] = (unsigned int *) calloc(d, sizeof(unsigned int));

        if(edges[i] == NULL)
            exit(EXIT_FAILURE);
    }

    x->label = -1;      // actually, since x->label is unsigned, this line performs the following assignment: x->label = ULLONG_MAX (however, that does not cause any problems)
    x->V = nodes;
    x->E = edges;
    x->pathsLength = 0;
}

void addGraph(graph* x) {
    unsigned int i, j, k;

    x->label++;

    char str[10 + 1], c;    // 2^32 - 1 = 4294967295 has ten digits

    for(i = 0; i < d; i++) {
        initializeVertex(&(x->V[i]));

        j = 0;

        while(j != d) {
            k = 0;

            c = getchar();
            while(c != ',' && c != '\n') {
                str[k] = c;
                k++;

                c = getchar();
            }

            str[k] = '\0';

            x->E[i][j] = stringToNumber(str);
            j++;
        }
    }

    x->pathsLength = 0;
}

void freeGraph(graph* x) {
    free(x->V);
    x->V = NULL;

    for(unsigned int i = 0; i < d; i++) {
        free(x->E[i]);
        x->E[i] = NULL;
    }
    free(x->E);
    x->E = NULL;
}

unsigned int parent(unsigned int index) {
    return index / 2;
}

unsigned int left(unsigned int index) {
    return 2 * index + 1;
}

unsigned int right(unsigned int index) {
    return 2 * index + 2;
}

void initializeMaxHeap(heap* x) {
    cell* A = (cell *) calloc(k, sizeof(cell));

    if(A == NULL)
        exit(EXIT_FAILURE);

    x->heapSize = 0;
    x->array = A;
}

void swapCells(cell* y, cell* z) {
    cell tmp;

    tmp = *y;
    *y = *z;
    *z = tmp;
}

void maxHeapify(heap* x, unsigned int index) {
    unsigned int l = left(index), r = right(index), max = index;

    if(l <= x->heapSize && x->array[l].length >= x->array[index].length)    // the use of greater than or equal sign is justified by the precedence rule (since this function is always called after extracting the current heap maximum)
        max = l;

    if(r <= x->heapSize && x->array[r].length >= x->array[max].length)      // the reasoning behind this logical expression can be seen above
        max = r;

    if(max != index) {
        swapCells(&(x->array[index]), &(x->array[max]));
        maxHeapify(x, max);
    }
}

void heapExtractMax(heap* x) {
    if(x->heapSize < 1)
        exit(EXIT_FAILURE);

    x->array[0] = x->array[x->heapSize - 1];    // set the last element of the array as the first one, basically overwriting the latter (i.e., the current maximum)
    x->heapSize--;      // decrease the overall size as one cell has been invalidated
    maxHeapify(x, 0);
}

void maxHeapRaiseUp(heap* x, unsigned int index) {
    while(index > 0 && x->array[parent(index)].length <= x->array[index].length) {      // it is requested that older graphs with same metric value as the newest added one should take precedence in the leaderboard (i.e., they must have "less chance" to be deleted with respect to the latest inserted graph)
        swapCells(&(x->array[index]), &(x->array[parent(index)]));
        index = parent(index);
    }
}

void maxHeapInsert(heap* x, cell n) {
    x->heapSize++;      // create a valid cell
    x->array[x->heapSize - 1] = n;      // set new element
    maxHeapRaiseUp(x, x->heapSize - 1);
}

void freeHeap(heap* x) {
    free(x->array);
    x->array = NULL;
}

void updateLeaderboard(heap* h, graph* g) {
    if(h->heapSize == k) {
        if(h->array[0].length > g->pathsLength)
            heapExtractMax(h);
    }

    if(h->heapSize < k) {
        cell newEntry;

        newEntry.tag = g->label;
        newEntry.length = g->pathsLength;

        maxHeapInsert(h, newEntry);
    }
}

void initializeQueue(queue* x) {
    unsigned int* indexes = (unsigned int *) calloc(d, sizeof(unsigned int));

    if(indexes == NULL)
        exit(EXIT_FAILURE);

    unsigned int* positions = (unsigned int *) calloc(d, sizeof(unsigned int));

    if(positions == NULL)
        exit(EXIT_FAILURE);

    x->heapSize = 0;
    x->graphVertices = NULL;
    x->vertexIndex = indexes;
    x->indexPositions = positions;
}

void overwriteQueue(queue* q, graph* g) {
    q->heapSize = 0;
    q->graphVertices = g->V;
}

void swap(unsigned int* y, unsigned int* z) {
    unsigned int tmp;

    tmp = *y;
    *y = *z;
    *z = tmp;
}

void minHeapify(queue* x, unsigned int index) {
    unsigned int l = left(index), r = right(index), min = index;

    if(l <= x->heapSize && x->graphVertices[x->vertexIndex[l]].distance <= x->graphVertices[x->vertexIndex[index]].distance)
        min = l;

    if(r <= x->heapSize && x->graphVertices[x->vertexIndex[r]].distance <= x->graphVertices[x->vertexIndex[min]].distance)
        min = r;

    if(min != index) {
        swap(&(x->vertexIndex[index]), &(x->vertexIndex[min]));
        swap(&(x->indexPositions[x->vertexIndex[index]]), &(x->indexPositions[x->vertexIndex[min]]));
        minHeapify(x, min);
    }
}

unsigned int heapExtractMin(queue* x) {
    if(x->heapSize < 1)
        exit(EXIT_FAILURE);

    unsigned int min = x->vertexIndex[0];

    x->vertexIndex[0] = x->vertexIndex[x->heapSize - 1];    // set the last element of the array as the first one, basically overwriting the latter (i.e., the current minimum)
    x->indexPositions[x->heapSize - 1] = 0;     // update the new first element's index
    x->heapSize--;      // invalidate previous minimum
    minHeapify(x, 0);

    return min;
}

void minHeapInsert(queue* x, unsigned int key) {
    x->heapSize++;
    x->vertexIndex[x->heapSize - 1] = key;
}

void minHeapUpdateKey(queue* x, unsigned int nodeIndex) {
    unsigned int nodePosition = x->indexPositions[nodeIndex];

    while(nodePosition > 0 && x->graphVertices[x->vertexIndex[parent(nodePosition)]].distance > x->graphVertices[x->vertexIndex[nodePosition]].distance) {
        swap(&(x->vertexIndex[nodePosition]), &(x->vertexIndex[parent(nodePosition)]));
        swap(&(x->indexPositions[nodeIndex]), &(x->indexPositions[x->vertexIndex[nodePosition]]));      // notice that, after the first swap, x->vertexIndex[parent(nodePosition)] has been changed into x->vertexIndex[nodePosition]
        nodePosition = parent(nodePosition);
    }
}

void freeQueue(queue* x) {
    free(x->vertexIndex);
    x->vertexIndex = NULL;

    free(x->indexPositions);
    x->indexPositions = NULL;
}

void optimizedDijkstra(graph* g, queue* q) {
    unsigned int i;

    overwriteQueue(q, g);

    g->V[0].distance = 0;
    g->V[0].isConnectedToSource = 1;

    for(i = 0; i < d; i++) {
        q->indexPositions[i] = i;
        minHeapInsert(q, i);
    }

    unsigned int minIndex;
    unsigned long long int neighbourVertexDistance;

    while(q->heapSize != 0) {
        minIndex = heapExtractMin(q);   // continue extracting until the queue is empty

        for(i = 0; i < d; i++) {    // i: used as the column entry of the adjacency matrix
            if(g->E[minIndex][i] != 0) {    // check only the adjacent nodes
                if(g->V[minIndex].isConnectedToSource == 1)
                    g->V[i].isConnectedToSource = 1;
                // neighbourVertexDistance = g->V[minIndex].distance + g->E[minIndex][i];
                neighbourVertexDistance = g->E[minIndex][i];
                if(g->V[minIndex].distance != ULLONG_MAX && g->V[minIndex].distance != 0)
                    neighbourVertexDistance += g->V[minIndex].distance;     // g->V[minIndex] := u, g->E[minIndex][i] := ArcWeight(u, v)

                if(g->V[i].distance > neighbourVertexDistance) {    // g->V[i] := v
                    g->V[i].distance = neighbourVertexDistance;
                    g->V[i].previous = &(g->V[minIndex]);
                    minHeapUpdateKey(q, i);     // update queue structure accordingly
                }
            }
        }

        if(g->V[minIndex].distance != ULLONG_MAX && g->V[minIndex].isConnectedToSource == 1)
            g->pathsLength += g->V[minIndex].distance;
    }
}

void showBestGraphs(heap x) {
    unsigned int i;

    if(x.heapSize != 0) {
        for(i = 0; i < k - 1 && i < x.heapSize - 1; i++)
            printf("%llu ", x.array[i].tag);
        printf("%llu\n", x.array[i].tag);
    } else
        printf("\n");
}

void parseInput() {
    if(scanf("%u %u\n", &d, &k) != 2)
        exit(EXIT_FAILURE);

    graph G;
    initializeGraph(&G);

    heap H;
    initializeMaxHeap(&H);

    queue Q;
    initializeQueue(&Q);

    char command[13 + 1 + 1];   // the longest input string, which is "AggiungiGrafo\n\0", has 15 characters
    short int code = !(EXIT_FAILURE);   // exit code

    if(fgets(command, sizeof(command), stdin) == NULL)
        code = EXIT_FAILURE;

    while(!feof(stdin) && code != EXIT_FAILURE) {
        if(command[0] == 'A') {   // insert new graph
            addGraph(&G);
            optimizedDijkstra(&G, &Q);
            updateLeaderboard(&H, &G);
        } else if (command[0] == 'T')     // print current leaderboard
            showBestGraphs(H);

        if(fgets(command, sizeof(command), stdin) == NULL)
            code = EXIT_FAILURE;
    }

    freeGraph(&G);
    freeHeap(&H);
    freeQueue(&Q);
}

int main() {
    parseInput();

    return 0;
}
